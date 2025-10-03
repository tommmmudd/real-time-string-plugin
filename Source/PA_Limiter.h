#pragma once

#define pa5_dlimlength 20

// ------------------------------------------------------------------------------------------
// A limiter, with fixed delayline.
// ------------------------------------------------------------------------------------------

/// Craig Webb limiter
class PA_Limiter
{
    
public:
    
    PA_Limiter() { }
    
    
    void reset()
    {
        for ( int i = 0; i < pa5_dlimlength; ++i )
        {
            m_dL[i] = 0.0;
            m_dR[i] = 0.0;
        }
        
        m_ptr   = 0;
        m_xpeak = 0.0;
        m_g     = 1.0;
        
    }
    
    // This returns the amount of gain reduction, so I can use it in the UI.
    float update(double inputL, double inputR, double& outputL, double& outputR)
    {
        double aL = std::fabs(inputL);
        double aR = std::fabs(inputR);
        double a  = std::fmax(aL, aR);
        
        double coeff;
        
        if ( a > m_xpeak )
            coeff = m_at;
        else
            coeff = m_rt;
        
        m_xpeak = (1.0 - coeff) * m_xpeak + coeff * a;
        
        double f = std::fmin(1.0, (m_lt / m_xpeak) );
        
        if (f < m_g)
            coeff = m_at;
        else
            coeff = m_rt;
        
        m_g = (1.0-coeff) * m_g + coeff * f;
        
        double dxL = m_dL[m_ptr];
        m_dL[m_ptr] = inputL;
        
        double dxR = m_dR[m_ptr];
        m_dR[m_ptr] = inputR;
        
        ++m_ptr;
        
        if ( m_ptr == pa5_dlimlength  )
        {
            m_ptr = 0;
        }
        
        outputL = dxL * m_g;
        outputR = dxR * m_g;
        
        float fmg = (float)m_g;
        
        return fmg;
    }
    
    
private:
    
    double m_dL[pa5_dlimlength];
    double m_dR[pa5_dlimlength];
    double m_xpeak;
    double m_g;
    double m_lt = 0.9;
    double m_at = 0.1;
    double m_rt = 0.0001;
    int   m_ptr;
    
};


