/*
  ==============================================================================

    GuitarString.h
    Created: 26 Jul 2024 10:50:29am
    Author:  Tom Mudd

  ==============================================================================
*/

#pragma once
#include <JuceHeader.h>
#include "PA_Limiter.h"


/// replace me!
class LowpassLazy {
    
public:
    
    void setSampleRate(float _sampleRate)
    {
        sampleRate = _sampleRate;
    }
    
    /// cutoff in Hz -
    void setCutoff(float newCutoff)
    {
        cutoff = newCutoff;
    }
    
    /// return lowpassed value. Call setSampleRate and setCutoff first
    float process(float newVal)
    {
        float cutoff_ratio = cutoff / sampleRate;
        float outVal = previousVal + cutoff_ratio * (newVal - previousVal);
        previousVal = outVal;  // update the previousVal with the new value
        return outVal;
    }
    
    
private:
    float previousVal = 0;
    float cutoff;
    float sampleRate = 44100;
    
};

/// Feedback with delay line, applying input to an interpolated point on the string
class FeedbackElement
{
public:
    FeedbackElement()
    {
        feedbackVector.resize(8000);
        lowpass.setCutoff(3500.0);
    }
    /// set position for feedback to act back on string (0-1)
    /// NEED TO HAVE ALREADY CALLED initialiseStringWithCurrentSettings()
    void setFeedbackPosition(float fbPos);
    
    void setFilterCutoff(float _newCutoff) { lowpass.setCutoff(_newCutoff); }
    
    /// call once size is initialised, before renderNextBlock
    void setNandH(int _N, float _h) { N = _N; h = _h; }
    
    void setFeedbackGain(float fbGain)  { feedbackGain = fbGain * 0.01; }
    
    void applyFeedback(double* u1, double* u2);
    
    void update(float newVal)
    {
        feedbackVector[fbWritePos] = lowpass.process(newVal * feedbackGain);
        fbWritePos++;
        fbWritePos %= feedbackVector.size();
    }
    
private:
    
    LowpassLazy lowpass;
    
    std::vector<float> feedbackVector;
    int N;
    float h;
    float feedbackGain = 0.001;
    float feedbackExcitationPoint = 0.8;
    int feedbackExcitationPointInterp;
    float feedbackExcitationPointFrac;
    int fbReadPos = 1;
    int fbWritePos = 0;
};






class GuitarString
{
public:
    GuitarString(int _whichString)
    {
        whichString = _whichString;
        i_s = std::to_string(whichString);
        fingerPosSmooth.setCurrentAndTargetValue(0.25);
        fingerForceSmooth.setCurrentAndTargetValue(0.0);
    }
    
    
    /// setup all derived params and other string settings for given N
    void setup();
    
    /// process a block of samples
    template <class FloatType>
    void processBlock(FloatType* outL, FloatType* outR, int numSamples, bool isStereo);

    /// rebuilds string (may be time consuming?
    void setSampleRate(double sampleRate)
    {
        SR = sampleRate;
        fingerPosSmooth.reset(SR, glideTime);
        fingerForceSmooth.reset(SR, 0.01);
        setup();
        limiter.reset();
    }
    
    /// position on string 0-1, or better 0.5-0.999
    void setFeedbackPosition(float _fbPos)
    {
        fbPos = _fbPos;
        
         if (isSetup)
         {
             feedbackElement.setFeedbackPosition(fbPos);
             externalFeedback.setFeedbackPosition(fbPos);
         }
         
    }
    
    /// feedback gain for the feedback element: e.g. 0-5
    void setInternalFeedbackGain(float _fbGain)
    {
        feedbackElement.setFeedbackGain(_fbGain);
    }
    
    /// feedback gain for the feedback element: e.g. 0-5
    void setExternalFeedbackGain(float _fbGain)
    {
        feedbackElement.setFeedbackGain(_fbGain);
    }
    
    
    void pluck(double _pluckPosition, double _pluckForce, double _pluckDuration=0.001)
    {
        isPlucking = true;
        pluckForce = _pluckForce;
        pluckPos = _pluckPosition;
        pluckCounter = 0;
        dur       = static_cast<int>(std::floor(_pluckDuration  * SR));
        ddur      = static_cast<double>(dur);

        xiint = (int)std::floor(pluckPos*L/h) - 1;
    }
    
    void inputExternalFeedback(const float* inputAudio, int numSamples)
    {
        for (int i=0; i<numSamples; i++)
        {
            externalFeedback.update(inputAudio[i]);
        }
    }
    

    
    /// link to APVTS
    void connectParams(juce::AudioProcessorValueTreeState& apvts)
    {
        fingerPosParam      = apvts.getRawParameterValue("fingerPos"+i_s);
        fingerForceParam    = apvts.getRawParameterValue("fingerForce"+i_s);
        outputPos1Param      = apvts.getRawParameterValue("outputPos1"+i_s);
        outputPos2Param      = apvts.getRawParameterValue("outputPos2"+i_s);
        
        externalFeedbackGainParam   = apvts.getRawParameterValue("externalFeedbackGain"+i_s);
        internalFeedbackGainParam   = apvts.getRawParameterValue("feedbackGain"+i_s);
        feedbackPosParam    = apvts.getRawParameterValue("feedbackPos"+i_s);
        feedbackFilterParam = apvts.getRawParameterValue("feedbackFilter"+i_s);
        
        stringLengthParam   = apvts.getRawParameterValue("stringLength"+i_s);
        stringTensionParam  = apvts.getRawParameterValue("stringTension"+i_s);
        stringRadiusParam   = apvts.getRawParameterValue("stringRadius"+i_s);
        stringT60_0Param    = apvts.getRawParameterValue("stringT60_0"+i_s);
        stringT60_1Param    = apvts.getRawParameterValue("stringT60_1"+i_s);
        
        fretHeightParam     = apvts.getRawParameterValue("fretHeight"+i_s);
        barrierHeightParam  = apvts.getRawParameterValue("barrierHeight"+i_s);
        
        barrierAlphaParam   = apvts.getRawParameterValue("barrierAlpha"+i_s);
        fretAlphaParam      = apvts.getRawParameterValue("fretAlpha"+i_s);
        fingerAlphaParam    = apvts.getRawParameterValue("fingerAlpha"+i_s);
        
        barrierStiffnessParam = apvts.getRawParameterValue("barrierStiffness"+i_s);
        fretStiffnessParam = apvts.getRawParameterValue("fretStiffness"+i_s);
        fingerStiffnessParam = apvts.getRawParameterValue("fingerStiffness"+i_s);
        
        fingerMassParam = apvts.getRawParameterValue("fingerMass"+i_s);
        
    }
    
    bool barrierOn    = true;
    bool kirchOn      = true;
    bool fretsOn      = true;
    bool fingerOn     = true;
    
    double* getString()         { return u; }
    int* getStringLength()      { return &N; }
    
    double* getFretPositions()  { return fretPositions; }
    int getNumFrets()           { return nFrets; }
    double* getFretHeight()     { return &fretsH; }
    double* getBarrierHeight()  { return &barrierH; }
    
    int status = 0; // 0 = active, 1 = explosion (zerod)
    
private:
    
    int whichString;    // is this string 1 or string 2?
    std::string i_s;
    
    // Explosion checking - trigger if output is above this value:
    double explosionThreshold = 200;
    PA_Limiter limiter;
    
    FeedbackElement feedbackElement;
    float fbPos = 0.5;
    
    FeedbackElement externalFeedback;
    
    bool isSetup    = false;
    
    // ==================================
    // finger parameters
    
    std::atomic<float>* fingerPosParam;
    std::atomic<float>* fingerForceParam;
    std::atomic<float>* outputPos1Param;
    std::atomic<float>* outputPos2Param;
    std::atomic<float>* externalFeedbackGainParam;
    std::atomic<float>* internalFeedbackGainParam;
    std::atomic<float>* feedbackPosParam;
    float currentFeedbackPos = 0.5;
    std::atomic<float>* feedbackFilterParam;
    float currentFeedbackFilter = 3500;
    
    
        
    std::atomic<float>* stringLengthParam;
    std::atomic<float>* stringTensionParam;
    std::atomic<float>* stringRadiusParam;
    std::atomic<float>* stringT60_0Param;
    std::atomic<float>* stringT60_1Param;
    
    std::atomic<float>* fretHeightParam;
    std::atomic<float>* barrierHeightParam;
    
    
    std::atomic<float>* barrierAlphaParam;
    std::atomic<float>* fretAlphaParam;
    std::atomic<float>* fingerAlphaParam;
    
    std::atomic<float>* barrierStiffnessParam;
    std::atomic<float>* fretStiffnessParam;
    std::atomic<float>* fingerStiffnessParam;
    
    std::atomic<float>* fingerMassParam;
    
    juce::SmoothedValue<double> fingerPosSmooth;
    juce::SmoothedValue<double> fingerForceSmooth;
    float glideTime = 0.1;
    
    double fingerPosValue = 0.25;
    double fingerForceValue = 0;
    
    // ===================================
    
    double pi         = juce::MathConstants<float>::pi;
    int Nf            = 44100*1;
    double SR         = 44100.0;
    
    
    
    double shiftV     = 2.220446049250313e-16;
    
    // barrier
    double barrierH   = - 0.001;                     // position of the barrier on the y axis in m
    double Kb         = 1.0e15;                      //  backboard stiffness
    double alphab     = 2.3;
    
    // frets
    const static int nFrets  = 20;                   // number of frets @TODO: unfreeze
    double STstretch  = 1.0;                         // fret stretch factor (semitones)
    double fretsH     = -0.0005;                     // fret height
    double Kf         = 1.0e15;                      // fret stiffness
    double alphaf     = 2.3;                         // fret alpha
    
    // finger
    double fForce     = 0.1;
    double fingPos1   = 0.07;
    double fingPos2   = 0.23;
    double w0         = 0.0;                         // initial position of the finger (y axis)
    double v0         = 0.0;                         // initial velocity
    double Mt         = 0.005;                       // finger mass (kg)
    double Kt         = 1.0e10;                      // finger collision stifness
    double alphat     = 1.3;                         // finger collision alpha
    
    // stability
    double epsilon    = 0.1;                         // Deviation from stability condition. epsilon > 0
    
    // Pluck
    double pluckForce = 0.2;
    double pluckTime  = 0.001;
    double pluckPos   = 0.81;
    double outposL    = 0.72;
    double outposR    = 0.82;
    
    // DC Blocker
    double xmL = 0.0;
    double ymL = 0.0;
    double xmR = 0.0;
    double ymR = 0.0;
    
    // ----------------------------------------------------------------------------------
    // Derived
    
    double rho    = 7850.0;
    double T0     = 12.1;    // 12.1,   12.3,    21.9,    39.2,    27.6,   49.2
    double radius = 0.00020; // 0.0002, 0.00015, 0.00015, 0.00015, 0.0001, 0.0001
    double E      = 2.0e11;
    double L      = 0.68;
    double T601   = 15.0;
    double T602   = 5.0;
    double Area;
    double rA;
    double I;
    double K;
    double c;
    double k;
    
    double *xt;
    double xtinc;
    double xtval;
    
    
    // ----------------------------------------------------------------------------------
    // Damping
    double omega1, omega2, param1, param2, sigma0, sigma1;
    
    // ----------------------------------------------------------------------------------
    // grid spacing
    double h;
    const static int maxN  = 300;
    int N  = 120;
    int nPoints = N-1;
    double lambda;
    
    int outpos1;
    int outpos2;
    
    // ----------------------------------------------------------------------------------
    // Memory
    alignas(32) double udata[maxN];
    alignas(32) double uPrevdata[maxN];
    alignas(32) double uNextdata[maxN];
    
    double *u     = udata;
    double *uPrev = uPrevdata;
    double *uNext = uNextdata;
    double *dummy_ptr;
    
    // for raised cosine plucking
    int dur;
    double ddur;
    int pluckCounter = 0;
    bool isPlucking = false;

    
    // double *out;
    
    int xiint;
    
    // ----------------------------------------------------------------------------------
    // Coeffs
    
    double mu, invM, tsk;
    
    double B0, BB0, B1, B2, C0, C1;
    
    double facBarr, facFrets, facFing;
    
    // ----------------------------------------------------------------------------------
    // Fretboard
    int ifr;
    
    double IF[nFrets * 2];
    int IFind[nFrets * 2];
    double fretPositions[nFrets * 2];
    
    int curr = 0;
    
    // ----------------------------------------------------------------------------------
    // More memory
    
    double wPrev, w, wNext;
    
    double pFac, fFac;

    alignas(32) double lin[maxN];
    alignas(32) double ITfull[maxN];
    alignas(32) double q[maxN];
    alignas(32) double b[maxN];
    alignas(32) double gb[maxN];
    alignas(32) double gf[maxN];
    alignas(32) double gt[maxN];
    alignas(32) double uDiff[maxN];
    alignas(32) double QG[maxN][4];
    alignas(32) double GQ[4][maxN];
    alignas(32) double As[4][4];
    alignas(32) double etaPosf[nFrets];
    alignas(32) double etaPosfAlpha[nFrets];
    alignas(32) double gradVf[nFrets];
    
    alignas(32) double gfTemp[maxN-1];
    alignas(32) double gfPrev[maxN-1];
    alignas(32) double gPrevb[maxN-1];
    alignas(32) double gbTemp[maxN-1];
    alignas(32) double gradVb[maxN-1];
    alignas(32) double uStr[maxN-1];
    alignas(32) double uStrPr[maxN-1];
    alignas(32) double linStr[maxN-1];
    alignas(32) double etaPosb[maxN-1];
    alignas(32) double etaPosbAlpha[maxN-1];
    
    
    double bs[4];
    double xs[4];
    
    double psiNextt, psiNextb, psiNextf, etaPost, psiPrevt, psiPrevb, psiPrevf;
    double facSAVfull, facSAVfullf, facSAV2full, facSAV2fullf, facKC;
    

    // ==============================================
    
    void updateString(double*u,  double* u1, double *u2, int N, double BB0, double B0, double B1,  double B2,  double C0, double C1)
    {
        int cp = 0;
        
        u[cp] = BB0*u1[cp] + B1*( u1[cp+1] ) + B2*( u1[cp+2] )
               + C0*u2[cp] + C1*( u2[cp+1] );
        
        ++cp;
        
        u[cp] = B0*u1[cp] + B1*( u1[cp-1] + u1[cp+1] ) + B2*( u1[cp+2] )
              + C0*u2[cp] + C1*( u2[cp-1] + u2[cp+1] );
        
        for( cp = 2; cp < N-3; ++cp)
        {
            u[cp] = B0*u1[cp] + B1*( u1[cp-1] + u1[cp+1] ) + B2*( u1[cp-2] + u1[cp+2] )
                  + C0*u2[cp] + C1*( u2[cp-1] + u2[cp+1] );
        }
        
        cp = N - 3;
        
        u[cp] = B0*u1[cp] + B1*( u1[cp-1] + u1[cp+1] ) + B2*( u1[cp-2] )
              + C0*u2[cp] + C1*( u2[cp-1] + u2[cp+1] );
        
        
        cp = N - 2;
        
        u[cp] = BB0*u1[cp] + B1*( u1[cp-1] ) + B2*( u1[cp-2] )
               + C0*u2[cp] + C1*( u2[cp-1] );
        
    }


    void solveGauss4(double A[4][4], double *b, double *x)
    {
        int s = 4;
        
        for (int j = 1; j < 4; ++j)
        {
            for (int i = s; i > j; --i)
            {
                double m = A[i-1][j-1] / A[j-1][j-1];
                
                for (int d = 0; d < s; ++d)
                {
                    A[i-1][d] -= m * A[j-1][d];
                }
                
                b[i-1] -= m*b[j-1];
            }
        }
        
        x[0] = 0.0;
        x[1] = 0.0;
        x[2] = 0.0;
        x[3] = b[3] / A[3][3];
        
        for (int i = s-1; i > 0; --i)
        {
            double sum = 0.0;
            
            for (int j = s; j > i; --j)
            {
                sum += A[i-1][j-1] * x[j-1];
            }
            
            x[i-1] = (b[i-1] - sum) / A[i-1][i-1];
        }
        
    }
    
};


