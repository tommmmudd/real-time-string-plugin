/*
  ==============================================================================

    DrawString.h
    Created: 27 Jul 2024 9:56:02am
    Author:  Tom Mudd

  ==============================================================================
*/

#pragma once

/**
 Visualise one or more GuitarStrings
 
 juce::Component for drawing the string profile, backboard and fretboard
 Currently set up to be able to draw two strings simultaneously
 */
class DrawString : public juce::Component, public juce::Timer
{
public:
    DrawString()
    {
        startTimer(40); // refresh rate
    }
    
    /// pass in a ptr to the string displacement array, and a pointer for the length
    void setPointerAndLength(double* _uPtr, int* _length)
    {
        ptrToString = _uPtr;
        stringLength = _length;
    }
    
    /// pass in a ptr to a second string displacement array, and a pointer for the length
    void setPointer2AndLength2(double* _uPtr, int* _length)
    {
        ptrToString2 = _uPtr;
        stringLength2 = _length;
    }
    
    /**
     Pass in a pointer to fret positions, along with info on the number of frets, and pointers to their height and the barrier height
     
     @param double* _fretPositions: a double array of fret positions
     @param int _numFrets: how many frets
     @param double* _fH ptr to fret height
     @param double* _bH ptr to backboard height
     */
    void setFretboard(double* _fretPositions, int _numFrets, double* _fH, double* _bH)
    {
        fretPositions = _fretPositions;
        numFrets = _numFrets;
        fretsH = _fH;
        barrierH = _bH;
    }
    
    void resized() override
    {
        
    }
    
    void paint(juce::Graphics& g) override
    {
        g.fillAll (juce::Colours::lightgrey);
        drawFrets(g);
        drawString(g);
    }
    
    void timerCallback() override
    {
        repaint();
    }
    
    /// draw the actual frets and the barrier
    void drawFrets(juce::Graphics& g)
    {
        if (numFrets > 0)
        {
            g.setColour (juce::Colours::purple);
            
            if (fretsOn)
            {
                for (int i=0; i<numFrets * 2; i+=2)
                {
                    // float fretPos = fretPositions[i];
                    // float strLen = float(*stringLength);
                    float pos = fretPositions[i] * getWidth(); //  / float(*stringLength);
                    juce::Line<float> line (juce::Point<float> (pos, scaleY(*fretsH)),  juce::Point<float> (pos, scaleY(*barrierH)));
                    g.drawLine (line, 2.0f);
                }
            }
            
            if (barrierOn)
            {
                juce::Line<float> line (juce::Point<float> (getWidth(), scaleY(*barrierH)),  juce::Point<float> (0, scaleY(*barrierH)));
                g.drawLine (line, 2.0f);
            }

            
        }
    }
    
    /// draw a single string profile (second string disabled in this version)
    void drawString(juce::Graphics& g)
    {
        g.setColour (juce::Colours::black);
        if (stringLength != nullptr)
        {
            // float pathY = getHeight() - 100;
            juce::Path path;
            
            path.startNewSubPath (juce::Point<float> (0, getHeight()*0.25));
            
            float dx = getWidth() / float(*stringLength);
            
            for (int i=0; i<(*stringLength)-1; i++)
            {
                double y = scaleY(ptrToString[i]);
                
                // NAN check
                if (y != y)
                    DBG("y point is NaN: " << y);
                else
                    path.lineTo (juce::Point<float> (dx * i, y));
            }
    
            //g.fillPath (path);
            g.strokePath (path, juce::PathStrokeType (3.0f));
        }
        
        
        /*
        g.setColour (juce::Colours::blue);
        if (stringLength2 != nullptr)
        {
            // float pathY = getHeight() - 100;
            juce::Path path;
            
            path.startNewSubPath (juce::Point<float> (0, getHeight()*0.25));
            
            float dx = getWidth() / float(*stringLength2);
            
            for (int i=0; i<(*stringLength2)-1; i++)
            {
                double y = scaleY(ptrToString2[i]);
                
                // NAN check
                if (y != y)
                    DBG("y point is NaN: " << y);
                else
                    path.lineTo (juce::Point<float> (dx * i, y));
            }
    
            //g.fillPath (path);
            g.strokePath (path, juce::PathStrokeType (3.0f));
        }
         */
    }
    
    float scaleY(float yZeroToOne)
    {
        return getHeight()*0.35 + yZeroToOne * -350 * getHeight();
    }
    
    bool fretsOn = true;
    bool barrierOn = true;
    
private:
    
    
    
    double* ptrToString = nullptr;
    double* ptrToString2 = nullptr;
    int* stringLength = nullptr;
    int* stringLength2 = nullptr;
    
    // int* fretboardIndexes;
    double* fretPositions;
    int numFrets = -1;
    double* fretsH = nullptr;
    double* barrierH = nullptr;
    
};
