/*
  ==============================================================================

    GuitarString.cpp
    Created: 26 Jul 2024 10:50:29am
    Author:  Tom Mudd

  ==============================================================================
*/

#include "GuitarString.h"

template <class FloatType>
void GuitarString::processBlock(FloatType* outL, FloatType* outR, int numSamples, bool isStereo)
{

    for (int n = 0; n < numSamples; ++n)
    {
        
        feedbackElement.setFeedbackGain(*internalFeedbackGainParam * 0.001);
        externalFeedback.setFeedbackGain(*externalFeedbackGainParam * 0.001);
        
        if (currentFeedbackPos != (*feedbackPosParam))
        {
            currentFeedbackPos = *feedbackPosParam;
            feedbackElement.setFeedbackPosition(currentFeedbackPos);
            externalFeedback.setFeedbackPosition(currentFeedbackPos);
        }
        
        if (currentFeedbackFilter != (*feedbackFilterParam))
        {
            currentFeedbackFilter = *feedbackFilterParam;
            feedbackElement.setFilterCutoff(currentFeedbackFilter);
            externalFeedback.setFilterCutoff(currentFeedbackFilter);
        }
        
        feedbackElement.applyFeedback(u, uPrev);
        externalFeedback.applyFeedback(u, uPrev);
        
        // Linear
        updateString(lin, u, uPrev, N, BB0, B0, B1, B2, C0, C1);
        lin[N-1] = 2.0 * u[N-1] - uPrev[N-1];
        
        // Add forces
        if (isPlucking)
        {
            int i_d = pluckCounter;
            double finSample = 0.5 * pluckForce * (1.0 - std::cos(M_PI * i_d / ddur));
            pluckCounter++;
            if (pluckCounter >= dur)
                isPlucking = false;

            lin[xiint] += pFac * finSample;         //xiint is pluck position as integer
        }
        
        if (fingerOn)
            lin[N-1]   += fFac * fingerForceValue; // fFac * fForce;
        
        
        // ======================================
        // Real time parameters
        
        // update finger mass (are others possible without resetting string?)
        Mt = *fingerMassParam;
        facSAVfullf  = 0.25*k*k/Mt;         // only two elements to update from finger mass?
        facSAV2fullf = k*k/Mt;
        
        // ALPHAS
        alphab = *barrierAlphaParam;
        alphaf = *fretAlphaParam;
        alphat = *fingerAlphaParam;
        
        // STIFFNESSES
        
        // barrier stiffness
        Kb = pow(10, *barrierStiffnessParam);
        facBarr  = h * Kb / (alphab + 1.0);
        
        // fret stiffness
        Kf = pow(10, *fretStiffnessParam);
        facFrets = Kf / (alphaf + 1.0);
        
        // finger stiffness
        Kt = pow(10, *fingerStiffnessParam);
        facFing  = Kt / (alphat + 1.0);
        
        
        
        // ======================================
        
        
        outposL   = *outputPos1Param;
        outposR   = *outputPos2Param;
        outpos1    = (int)std::floor(outposL * N);
        outpos2    = (int)std::floor(outposR * N);
        
        for (int i = 0; i < nPoints; ++i)
        {
            uStr[i]   = u[i];
            uStrPr[i] = uPrev[i];
            linStr[i] = lin[i];
        }
        
        
        if (kirchOn)
        {
            double ih = 1.0 / (h * h);
            
            q[0] = (-2.0 * uStr[0] + uStr[1]) * ih;
            
            for (int cp = 1; cp < nPoints-1; ++cp)
            {
                q[cp] = (-2.0 * uStr[cp] + uStr[cp+1] + uStr[cp-1]) * ih;
            }
            
            q[nPoints - 1] = (-2.0 * uStr[nPoints - 1] + uStr[nPoints - 2]) * ih;
            
            q[nPoints] = 0.0;
        }
        
        
        if (barrierOn)
        {
            for (int i = 0; i < nPoints; ++i)
            {
                double etab = -uStr[i] + barrierH;
                
                if (etab > 0.0)
                {
                    etaPosb[i] = etab;
                }
                else
                {
                    etaPosb[i] = 0.0;
                }
                
                etaPosbAlpha[i] = std::pow(etaPosb[i], alphab);
                
                gPrevb[i] = gb[i];
            }
            
            double sum = 0.0;
            
            for (int i = 0; i < nPoints; ++i)
            {
                sum += etaPosbAlpha[i] * etaPosb[i];
            }
            
            double Vb    = facBarr * sum;
            double invVb = 1.0/std::sqrt(2.0*(Vb + shiftV));
            
            for (int i = 0; i < nPoints; ++i)
            {
                gradVb[i] = -h*Kb*etaPosbAlpha[i];
                gbTemp[i] = invVb*gradVb[i];
            }
            
            double zetab = 0.0;
            
            for (int i = 0; i < nPoints; ++i)
            {
                zetab += gbTemp[i] * (linStr[i] - uStrPr[i]);
            }
            
            if (zetab < 0.0 || zetab > 0.0)
            {
                if (zetab < -4.0*psiPrevb)
                {
                    double gammab = -4.0*psiPrevb/zetab;
                    
                    for (int i = 0; i < nPoints; ++i)
                    {
                        gbTemp[i] = gammab*gbTemp[i];
                    }
                }
            }
            else
            {
                double zPrevb = 0.0;
                    
                for (int i = 0; i < nPoints; ++i)
                {
                    zPrevb += gPrevb[i] * (linStr[i] - uStrPr[i]);
                }
                    
                if (zPrevb < 0.0 || zPrevb > 0.0)
                {
                    double gammab = -4.0 * psiPrevb / zPrevb;
                    
                    if (gammab < 0.0)
                    {
                        gammab = 0.5*gammab;
                    }
                    else
                    {
                        gammab = 1.5*gammab;
                    }
                        
                    for (int i = 0; i < nPoints; ++i)
                    {
                        gbTemp[i] = gammab * gPrevb[i];
                    }
                }
            }
            
            
            for (int i = 0; i < nPoints; ++i)
            {
                gb[i] = gbTemp[i];
            }
            gb[N-1] = 0.0;
        }
        
        
        
        if (fretsOn)
        {
            for (int i = 0; i < nFrets; ++i)
            {
                double temp = uStr[IFind[i*2]-1] * IF[i*2];
                temp       += uStr[IFind[i*2+1]-1] * IF[i*2+1];
                double etaf = fretsH - temp;
                
                if (etaf > 0.0)
                {
                    etaPosf[i] = etaf;
                }
                else
                {
                    etaPosf[i] = 0.0;
                }
                
                etaPosfAlpha[i] = std::pow(etaPosf[i], alphaf);
            }
            
            for (int i = 0; i < nPoints; ++i)
            {
                gfPrev[i] = gf[i];
            }
            
            double sum = 0.0;
            
            for (int i = 0; i < nFrets; ++i)
            {
                sum += etaPosfAlpha[i]*etaPosf[i];
            }
            
            double Vf = facFrets * sum;
            double invVf = 1.0/std::sqrt(2.0*(Vf + shiftV));
            
            for (int i = 0; i < nFrets; ++i)
            {
                gradVf[i] = -Kf*etaPosfAlpha[i];
            }
            
            for (int i = 0; i < nPoints; ++i)
            {
                gfTemp[i] = 0.0;
            }
            
            for (int i = 0; i < nFrets; ++i)
            {
                int ind1  = IFind[i*2]-1;
                double v1 = IF[i*2];
                int ind2 = IFind[i*2+1]-1;
                double v2 = IF[i*2+1];
                
                gfTemp[ind1] += v1 * invVf * gradVf[i];
                gfTemp[ind2] += v2 * invVf * gradVf[i];
            }
            
            double zetaf = 0.0;
            
            for (int i = 0; i < nPoints; ++i)
            {
                zetaf += gfTemp[i] * (linStr[i] - uStrPr[i]);
            }
            
            if (zetaf < 0.0 || zetaf > 0.0)
            {
                if (zetaf < -4.0*psiPrevf)
                {
                    double gammaf = -4.0*psiPrevf/zetaf;
                    
                    for (int i = 0; i < nPoints; ++i)
                    {
                        gfTemp[i] = gammaf * gfTemp[i];
                    }
                }
            }
            else
            {
                double zPrevf = 0.0;
                
                for (int i = 0; i < nPoints; ++i)
                {
                    zPrevf += gfPrev[i] * (linStr[i] - uStrPr[i]);
                }
                
                if (zPrevf < 0.0 || zPrevf > 0.0)
                {
                    double gammaf = -4.0 * psiPrevf / zPrevf;
                    
                    if (gammaf < 0.0)
                    {
                        gammaf = 0.5*gammaf;
                    }
                    else
                    {
                        gammaf = 1.5*gammaf;
                    }
                    
                    for (int i = 0; i < nPoints; ++i)
                    {
                        gfTemp[i] = gammaf * gfPrev[i];
                    }
                }
                
            }
            
            for (int i = 0; i < nPoints; ++i)
            {
                gf[i] = gfTemp[i];
            }
            gf[N-1] = 0.0;
        }
        
        // @TODO: add finger placement in xt[] externally
        
        
        // update finger position even if not on string
        fingerPosSmooth.setTargetValue(*fingerPosParam);
        fingerForceSmooth.setTargetValue(*fingerForceParam * *fingerForceParam * 3.0);
        
        fingerPosValue = fingerPosSmooth.getNextValue();
        fingerForceValue = fingerForceSmooth.getNextValue();
        
        
        
        if (fingerOn)
        {
            // scaled up to 0-5 cubic
            fingerForceValue = (fingerForceValue * fingerForceValue * fingerForceValue) * 5.0;
            
            for (int i = 0; i < N; ++i)
            {
                ITfull[i] = 0.0;
            }
            
            ITfull[N-1] = -1.0;
            
            
            // set position
            double indt = fingerPosValue * L / h;
            int xind = (int)std::floor(indt);
            double fract = indt - (double)xind;
            
            if (xind < 1)
            {
                xind = 1;
                fract = 0.0;
            }
            else if (xind > nPoints - 3)
            {
                xind = nPoints - 3;
                fract = 0.0;
            }
            
            ITfull[xind-2] = -fract*(fract-1.0)*(fract-2.0)/(6.0);
            ITfull[xind-1] = (fract-1.0)*(fract+1.0)*(fract-2.0)/2.0;
            ITfull[xind]   = -fract*(fract+1.0)*(fract-2.0)/2.0;
            ITfull[xind+1] = fract*(fract+1.0)*(fract-1.0)/6.0;
            
            double etat = 0.0;
            
            for (int i = 0; i < N; ++i)
            {
                etat += ITfull[i] * u[i];
            }
            
            etaPost = 0.5*(etat + std::abs(etat));
            
            double etaPostAlpha = std::pow(etaPost,alphat);
            double Vt = facFing*etaPostAlpha*etaPost;
            double invVt = 1.0/std::sqrt(2.0*(Vt + shiftV));
            double gradVt = Kt*etaPostAlpha;
            
            for (int i = 0; i < N; ++i)
            {
                gt[i] = ITfull[i] * gradVt * invVt;
            }
        }
        
        
        double b1 = 0.0;
        
        for (int i = 0; i < N; ++i)
        {
            b1 += gb[i] * uPrev[i];
        }
        
        
        double b2 = 0.0;
        
        for (int i = 0; i < N; ++i)
        {
            b2 += q[i] * uPrev[i];
        }
        
        
        double b3 = 0.0;
        
        for (int i = 0; i < N; ++i)
        {
            b3 += gf[i] * uPrev[i];
        }
        
        double b4 = 0.0;
        
        for (int i = 0; i < N; ++i)
        {
            b4 += gt[i] * uPrev[i];
        }
        
        
        for (int i = 0; i < nPoints; ++i)
        {
            double b_1 = facSAVfull*gb[i]*b1;
            double b_2 = facSAVfull*gf[i]*b3;
            double b_3 = facSAVfull*gt[i]*b4;
            double b_4 = facSAV2full*gb[i]*psiPrevb;
            double b_5 = facSAV2full*gf[i]*psiPrevf;
            double b_6 = facSAV2full*gt[i]*psiPrevt;
            double b_7 = facKC*q[i]*b2;
            
            b[i] = lin[i] + b_1 + b_2 + b_3 - b_4 - b_5 - b_6 - b_7;
            
        }
            
        b[N-1] = lin[N-1] + facSAVfullf*gt[N-1]*b4 - facSAV2fullf*gt[N-1]*psiPrevt;
    
        
        for (int i = 0; i < N-1; ++i)
        {
            QG[i][0] = facSAVfull*gb[i];
            QG[i][1] = facKC*q[i];
            QG[i][2] = facSAVfull*gf[i];
            QG[i][3] = facSAVfull*gt[i];
            GQ[0][i] = gb[i];
            GQ[1][i] = q[i];
            GQ[2][i] = gf[i];
            GQ[3][i] = gt[i];
        }
        
        QG[N-1][0] = 0.0;
        QG[N-1][1] = 0.0;
        QG[N-1][2] = 0.0;
        QG[N-1][3] = facSAVfullf*gt[N-1];
        GQ[0][N-1] = 0.0;
        GQ[1][N-1] = 0.0;
        GQ[2][N-1] = 0.0;
        GQ[3][N-1] = gt[N-1];
        
        
        for (int i = 0; i < 4; ++i)
        {
            for (int j = 0; j < 4; ++j)
            {
                double sum = 0.0;
                
                for (int k = 0; k < N; ++k)
                {
                    sum += GQ[i][k] * QG[k][j];
                }
                
                As[i][j] = sum;
            }
        }
        
        As[0][0] += 1.0;
        As[1][1] += 1.0;
        As[2][2] += 1.0;
        As[3][3] += 1.0;
        
        for (int i = 0; i < 4; ++i)
        {
            double sum = 0.0;
            
            for (int j = 0; j < N; ++j)
            {
                sum += GQ[i][j] * b[j];
            }
            bs[i] = sum;
        }
        
        solveGauss4(As, bs, xs);
        
        for (int i = 0; i < N; ++i)
        {
            uNext[i] = b[i] - (QG[i][0]*xs[0] + QG[i][1]*xs[1] + QG[i][2]*xs[2] + QG[i][3]*xs[3]);
        }
        
        for (int i = 0; i < N; ++i)
        {
            uDiff[i] = uNext[i] - uPrev[i];
        }
        
        
        double sumb = 0.0;
        
        for (int i = 0; i < N; ++i)
        {
            sumb += 0.5 * gb[i] * uDiff[i];
        }
        
        psiNextb = psiPrevb + sumb;

        sumb = 0.0;
        
        for (int i = 0; i < N; ++i)
        {
             sumb += 0.5 * gf[i] * uDiff[i];
        }
        
        psiNextf = psiPrevf + sumb;
        
        sumb = 0.0;
        
        for (int i = 0; i < N; ++i)
        {
            sumb += 0.5 * gt[i] * uDiff[i];
        }
        
        psiNextt = psiPrevt + sumb;
        
        psiPrevb = psiNextb;
        psiPrevf = psiNextf;
        psiPrevt = psiNextt;
        
        
        // ================================
        // explosion check
        
        // has actually achieved NaN?
        if (uNext[outpos1] != uNext[outpos1] || uNext[outpos2] != uNext[outpos2])
        {
            DBG("uNext[outpos1] is NaN " << uNext[outpos1] );
        }
        
        // Getting too big?
        else if (uNext[outpos1] * 500.0 > explosionThreshold || uNext[outpos2] * 500.0 > explosionThreshold)
        {
            DBG(" EXPLOSION " << uNext[outpos1] * 500.0);
            // outL[n] += 0;
            status = 1;
            
            feedbackElement.update(0);
        }
        
        // is already exploded?
        else if (status == 1)
        {
            // outL[n] += 0;
            feedbackElement.update(0);
        }
        
        // otherwise output as normal
        else
        {
            double outValLPreLimiter = uNext[outpos1] * 200.0;
            double outValRPreLimiter = uNext[outpos2] * 200.0;
            double outValLPostLimiter;
            double outValRPostLimiter;
            
            // Limiter
            limiter.update (outValLPreLimiter, outValRPreLimiter, outValLPostLimiter, outValRPostLimiter);
            
            auto outValL = static_cast<FloatType> (outValLPreLimiter);
            auto outValR = static_cast<FloatType> (outValRPreLimiter);
            
            outL[n] += outValL;
            if (isStereo)
                outR[n] += outValR;
            else
                outL[n] += outValR;
            
            feedbackElement.update ((outValL + outValR) * 0.5);
        }
        
        // swap pointers
        dummy_ptr = uPrev;
        uPrev = u;
        u     = uNext;
        uNext = dummy_ptr;
        
        
    }
    
    // DC Block
    for (int n = 0; n < numSamples; ++n)
    {
        auto xL = outL[n];
        auto yL = xL - xmL + 0.995 * ymL;

        outL[n] += yL;

        xmL = xL;
        ymL = yL;
        
        auto xR = outR[n];
        auto yR = xR - xmR + 0.995 * ymR;

        outR[n] += yR;

        xmR = xR;
        ymR = yR;
    }
    
    
    
}

void GuitarString::setup()
{
    
    // reenable output
    status = 0;
    
    // ==============================
    // updates from parameters
    
    // string params
    L = *stringLengthParam;
    T0 = *stringTensionParam;
    radius = *stringRadiusParam;
    T601 = *stringT60_0Param;
    T602 = *stringT60_1Param;
    
    // make sure high decay is always less than low decay
    if (T602 >= T601)
    {
        T602 = T601 - 1;
        DBG("High decay is higher then low decay. Setting to " << T602);
    }

    // fret and barrier
    fretsH = *fretHeightParam;
    barrierH = *barrierHeightParam;
    
    // alphas
    alphab = *barrierAlphaParam;
    alphaf = *fretAlphaParam;
    alphat = *fingerAlphaParam;
    
    // stiffness
    Kb = pow(10, *barrierStiffnessParam);
    Kf = pow(10, *fretStiffnessParam);
    Kt = pow(10, *fingerStiffnessParam);
    
    // finger mass
    Mt = *fingerMassParam;
    
    
    // ==============================
    Area   = pi*radius*radius;
    rA     = rho*Area;
    I      = (pi*radius*radius*radius*radius) / 4.0;
    K      = std::sqrt( E * I / rA );
    c      = std::sqrt( T0 / rA );
    k      = 1.0 / SR;
    
     
    
    // ----------------------------------------------------------------------------------
    // Damping
    omega1 = 0.0 * 2.0 * pi;
    omega2 = 1000.0 * 2.0 * pi;
    param1 = (-c*c + std::sqrt(c*c*c*c + 4.0 * K*K * omega1*omega1)) / (2.0 * K*K);
    param2 = (-c*c + std::sqrt(c*c*c*c + 4.0 * K*K * omega2*omega2)) / (2.0 * K*K);

    sigma0 = (6.0*std::log(10.0) / (param2 - param1))*(param2/T601 - param1/T602);
    sigma1 = (6.0*std::log(10.0) / (param2 - param1))*(-1.0/T601 + 1.0/T602);
    
    
    // ----------------------------------------------------------------------------------
    // Compute grid spacing
    h      = sqrt((c*c*k*k + 4.0*sigma1*k + std::sqrt(  ((c*c*k*k + 4.0*sigma1*k)
                    * (c*c*k*k + 4.0*sigma1*k)) + ((16.0*K*K*k*k))))/2.0) * (1.0 + epsilon);
    N         = (int)(std::floor(L/h));
    if (N > maxN)
        N = maxN;
    h             = L / ((double)N);
    nPoints   = N - 1;           // 119
    lambda = c*k/h;
    
    outpos1    = (int)std::floor(outposL * N);
    outpos2    = (int)std::floor(outposR * N);
    
    DBG("N = " << N);
    
    
    // ----------------------------------------------------------------------------------
    // Coeffs
    
    mu   = K * k / (h * h);
    invM = 1.0 / (1.0 + sigma0 * k);
    tsk  = (2.0 * sigma1 * k) / (h * h);
    
    B0   = invM * (2.0 - 2.0*lambda*lambda - 6.0*mu*mu - 2.0*tsk);
    BB0  = invM * (2.0 - 2.0*lambda*lambda - 5.0*mu*mu - 2.0*tsk);
    B1   = invM * (lambda*lambda + 4.0*mu*mu + tsk);
    B2   = invM * (- mu*mu);
    C0   = invM * ((-1.0 + sigma0*k) + 2.0*tsk);
    C1   = invM * -tsk;
    
    
    facBarr  = h * Kb / (alphab + 1.0);
    facFrets = Kf / (alphaf + 1.0);
    facFing  = Kt / (alphat + 1.0);
    
    // ----------------------------------------------------------------------------------
    // Fretboard
    
    // @TODO: unfix number of frets (dynamic arrays)
    // ifr = nFrets * 2;
    
    // IF[ifr];
    // IFind[ifr];
    
    curr = 0;
    
    for (int i = 0; i < nFrets; ++i)
    {
        int ip     = i + 1;
        double id  = (double)ip;
        double xST = L * (1.0 - std::pow(2.0, -id * STstretch / 12.0));
        
        int xind = (int)std::floor(xST / h);
        
        double fracf = xST / h - (double)xind;
        
        IF[curr]    = 1.0 - fracf;
        IFind[curr] = xind;
        
        // Tom: added this to have a single float array representing the (unsnapped) fret positions 0-1
        fretPositions[curr] = xST / L;
        
        ++curr;
        
        IF[curr]    = fracf;
        IFind[curr] = xind + 1;
        
        ++curr;
        
        
        // ----------------------------------------------------------------------------------
        // More memory
        
        
        wPrev = w0;
        w     = w0 + k*v0;
        wNext = 0.0;
        
        pFac  = (k*k/rA) / (1+sigma0*k) / h;
        fFac  = -k*k/Mt;                        // timestep^2 / finger mass (Mt)
        
        u[N-1] = w;
        uPrev[N-1] = wPrev;
        uNext[N-1] = wNext;
        
        for (int i = 0; i < N; ++i)
        {
            u[i]      = 0.0;
            uNext[i]  = 0.0;
            uPrev[i]  = 0.0;
            
            lin[i]    = 0.0;
            ITfull[i] = 0.0;
            q[i]      = 0.0;
            b[i]      = 0.0;
            gb[i]     = 0.0;
            gf[i]     = 0.0;
            gt[i]     = 0.0;
            uDiff[i]  = 0.0;
            
            QG[i][0] = 0.0;
            QG[i][1] = 0.0;
            QG[i][2] = 0.0;
            QG[i][3] = 0.0;
            GQ[0][i] = 0.0;
            GQ[1][i] = 0.0;
            GQ[2][i] = 0.0;
            GQ[3][i] = 0.0;
        }
        
        for (int i = 0; i < nPoints; ++i)
        {
            gfTemp[i] = 0.0;
            gfPrev[i] = 0.0;
            gPrevb[i] = 0.0;
            gbTemp[i] = 0.0;
            gradVb[i] = 0.0;
            uStr[i]   = 0.0;
            uStrPr[i] = 0.0;
            linStr[i] = 0.0;
            etaPosb[i]      = 0.0;
            etaPosbAlpha[i] = 0.0;
        }
        
        for (int i = 0; i < 4; ++i)
        {
            bs[i] = 0.0;
            xs[i] = 0.0;
        }
        
        psiNextt = 0.0;
        psiNextb = 0.0;
        psiNextf = 0.0;
        etaPost  = 0.0;
        psiPrevt = 0.0;
        psiPrevb = std::sqrt(2.0*(shiftV));
        psiPrevf = std::sqrt(2.0*(shiftV));
        
        facSAVfull   = invM * 0.25*k*k/(rA*h);
        facSAVfullf  = 0.25*k*k/Mt;
        facSAV2full  = invM * (k*k/(rA*h));
        facSAV2fullf = k*k/Mt;
        facKC = invM * k*k*(E*h*0.25/(L*rho));
    }
    
    isSetup = true;
    
    feedbackElement.setNandH(N, h);
    feedbackElement.setFeedbackPosition(currentFeedbackPos);
    externalFeedback.setNandH(N, h);
    externalFeedback.setFeedbackPosition(currentFeedbackPos);
}


template void GuitarString::processBlock(double* outL, double* outR, int numSamples, bool isStereo);
template void GuitarString::processBlock(float* outL,  float* outR,  int numSamples, bool isStereo);






void FeedbackElement::setFeedbackPosition(float fbPos)
{
    feedbackExcitationPoint = fbPos;
    feedbackExcitationPointInterp = 1 + floor(N * feedbackExcitationPoint);
    feedbackExcitationPointFrac = 1 + feedbackExcitationPoint / h - feedbackExcitationPointInterp;
}


void FeedbackElement::applyFeedback(double* u1, double* u2)
{
    
    // read from feedback vector and excite string at that point (better to do raised cosine?)
    u2[feedbackExcitationPointInterp] += (1 - feedbackExcitationPointFrac) * feedbackVector[fbReadPos];
    u2[feedbackExcitationPointInterp + 1] += feedbackExcitationPointFrac * feedbackVector[fbReadPos];
    u1[feedbackExcitationPointInterp] += (1 - feedbackExcitationPointFrac) * feedbackVector[fbReadPos];
    u1[feedbackExcitationPointInterp + 1] += feedbackExcitationPointFrac * feedbackVector[fbReadPos];
    fbReadPos ++;
    fbReadPos %= feedbackVector.size();
}
