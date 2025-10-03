/*
  ==============================================================================

    GuitarControlComponent.h
    Created: 11 Feb 2025 3:33:40pm
    Author:  Tom Mudd

  ==============================================================================
*/

#pragma once
#include <JuceHeader.h>
#include "DrawString.h"
#include "GuitarString.h"
#include "PluginProcessor.h"


class GuitarControlComponent : public juce::Component
{
public:
    
    GuitarControlComponent(StringRTAudioProcessor& p, GuitarString& gs, DrawString& ds, int ind);
    
    void paint(juce::Graphics& g) override;
    void resized() override;
    
    
    void setFlag(juce::Button* button, juce::String name)
    {
        auto state = button->getToggleState();
        if (name == "kirch")
            guitarString.kirchOn = state;
        else if (name == "barrier")
        {
            guitarString.barrierOn = state;
            drawString.barrierOn = state;
        }
        else if (name == "frets")
        {
            guitarString.fretsOn = state;
            drawString.fretsOn = state;
        }
        else if (name == "finger")
        {
            guitarString.fingerOn = state;
        }
    }
    
private:
    StringRTAudioProcessor& audioProcessor;
    GuitarString& guitarString;
    DrawString& drawString;
    
    int index;          // which string are we addressing? 1, 2, etc
    std::string i_s;    // that index as a string, "1", "2"
    
    juce::TextButton resetButton;
    juce::TextButton pluckButton;
    
    
    std::atomic<float>* pluckForceParam;
    std::atomic<float>* pluckPosParam;
    
    void putFingerAtFret(int newFret)
    {
        float positionOnFret = 0.5;   // 0.667 is 2/3 of the way towards higher fret
        
        int fretIndex = (newFret) * 2;
        auto fretPositions = guitarString.getFretPositions();
        
        // fret is 1-20
        float fingerPos = fretPositions[fretIndex] * positionOnFret;
        
        // if we're higher than first fret
        if (newFret > 1)
        {
            int previousIndex = fretIndex - 2;
            float nextFretPos = static_cast<float>(fretPositions[previousIndex]);
            float previousFretPos = static_cast<float>(fretPositions[previousIndex]);
            fingerPos = (nextFretPos - previousFretPos) * positionOnFret  +  previousFretPos;
        }
        
        fingerPosSlider.setValue(fingerPos);
    }
    
    
    // REAL TIME SLIDERS
    juce::Label fingerPosLabel { {}, "Finger Position (0-1)" };
    juce::Slider fingerPosSlider;
    std::unique_ptr<juce::AudioProcessorValueTreeState::SliderAttachment> fingerPosSliderAttachment;
    
    juce::Label fretPositionLabel { {}, "Finger at Fret" };
    juce::Slider fretPositionSlider;
    std::unique_ptr<juce::AudioProcessorValueTreeState::SliderAttachment> fretPositionSliderAttachment;

    
    juce::Label fingerForceLabel { {}, "Finger Force" };
    juce::Slider fingerForceSlider;
    std::unique_ptr<juce::AudioProcessorValueTreeState::SliderAttachment> fingerForceSliderAttachment;
    
    juce::Label outputPos1Label { {}, "Output Pos 1" };
    juce::Slider outputPos1Slider;
    std::unique_ptr<juce::AudioProcessorValueTreeState::SliderAttachment> outputPos1SliderAttachment;
    
    juce::Label outputPos2Label { {}, "Output Pos 2" };
    juce::Slider outputPos2Slider;
    std::unique_ptr<juce::AudioProcessorValueTreeState::SliderAttachment> outputPos2SliderAttachment;
    
    juce::Label externalFeedbackGainLabel { {}, "External Feedback Gain" };
    juce::Slider externalFeedbackGainSlider;
    std::unique_ptr<juce::AudioProcessorValueTreeState::SliderAttachment> externalFeedbackSliderAttachment;
    
    juce::Label internalFeedbackGainLabel { {}, "Internal Feedback Gain" };
    juce::Slider internalFeedbackGainSlider;
    std::unique_ptr<juce::AudioProcessorValueTreeState::SliderAttachment> internalFeedbackGainSliderAttachment;
    
    juce::Label feedbackPosLabel { {}, "Feedback Pos" };
    juce::Slider feedbackPosSlider;
    std::unique_ptr<juce::AudioProcessorValueTreeState::SliderAttachment> feedbackPosSliderAttachment;
    
    juce::Label feedbackFilterLabel { {}, "Feedback Filter" };
    juce::Slider feedbackFilterSlider;
    std::unique_ptr<juce::AudioProcessorValueTreeState::SliderAttachment> feedbackFilterSliderAttachment;
    
    juce::Label pluckForceLabel { {}, "Pluck Force" };
    juce::Slider pluckForceSlider;
    std::unique_ptr<juce::AudioProcessorValueTreeState::SliderAttachment> pluckForceSliderAttachment;
    
    juce::Label pluckPosLabel { {}, "Pluck Position" };
    juce::Slider pluckPosSlider;
    std::unique_ptr<juce::AudioProcessorValueTreeState::SliderAttachment> pluckPosSliderAttachment;
    
    juce::ToggleButton kirchButton   { "Kirchhoff" };
    juce::ToggleButton barrierButton   { "barrier" };
    juce::ToggleButton fretsButton   { "frets" };
    juce::ToggleButton fingerButton   { "finger" };
    
    std::unique_ptr<juce::AudioProcessorValueTreeState::ButtonAttachment> kcAttachment;
    std::unique_ptr<juce::AudioProcessorValueTreeState::ButtonAttachment> barrierAttachment;
    std::unique_ptr<juce::AudioProcessorValueTreeState::ButtonAttachment> fretsAttachment;
    std::unique_ptr<juce::AudioProcessorValueTreeState::ButtonAttachment> fingerAttachment;
    // juce::AudioProcessorValueTreeState::ButtonAttachment>
    
    juce::Label volumeLabel { {}, "Volume" };
    juce::Slider volumeSlider;
    std::unique_ptr<juce::AudioProcessorValueTreeState::SliderAttachment> volumeSliderAttachment;
    
    // status info label (e.g. explosion notifications)
    juce::Label statusLabel { {}, "Status:" };
    
    
    
    
    
    // STRING PARAMS
    juce::Label stringLengthLabel { {}, "String Length" };
    juce::Slider stringLengthSlider;
    std::unique_ptr<juce::AudioProcessorValueTreeState::SliderAttachment> stringLengthSliderAttachment;

    juce::Label stringTensionLabel { {}, "String Tension" };
    juce::Slider stringTensionSlider;
    std::unique_ptr<juce::AudioProcessorValueTreeState::SliderAttachment> stringTensionSliderAttachment;

    juce::Label stringRadiusLabel { {}, "String Radius" };
    juce::Slider stringRadiusSlider;
    std::unique_ptr<juce::AudioProcessorValueTreeState::SliderAttachment> stringRadiusSliderAttachment;
    
    juce::Label stringT60_0Label { {}, "String Low Decay" };
    juce::Slider stringT60_0Slider;
    std::unique_ptr<juce::AudioProcessorValueTreeState::SliderAttachment> stringT60_0SliderAttachment;
    
    juce::Label stringT60_1Label { {}, "String High Decay" };
    juce::Slider stringT60_1Slider;
    std::unique_ptr<juce::AudioProcessorValueTreeState::SliderAttachment> stringT60_1SliderAttachment;
    
    
    // FRET and Barrier PARAMS
    juce::Label fretHeightLabel { {}, "Fret height" };
    juce::Slider fretHeightSlider;
    std::unique_ptr<juce::AudioProcessorValueTreeState::SliderAttachment> fretHeightSliderAttachment;
    
    juce::Label barrierHeightLabel { {}, "Barrier height" };
    juce::Slider barrierHeightSlider;
    std::unique_ptr<juce::AudioProcessorValueTreeState::SliderAttachment> barrierHeightSliderAttachment;
    
    
    
    // ======================
    // alphas
    
    juce::Label alphaBarrierLabel { {}, "Barrier Alpha" };
    juce::Slider alphaBarrierSlider;
    std::unique_ptr<juce::AudioProcessorValueTreeState::SliderAttachment> alphaBarrierSliderAttachment;
    
    juce::Label alphaFretLabel { {}, "Fret Alpha" };
    juce::Slider alphaFretSlider;
    std::unique_ptr<juce::AudioProcessorValueTreeState::SliderAttachment> alphaFretSliderAttachment;
    
    juce::Label alphaFingerLabel { {}, "Finger Alpha (N)" };
    juce::Slider alphaFingerSlider;
    std::unique_ptr<juce::AudioProcessorValueTreeState::SliderAttachment> alphaFingerSliderAttachment;
        
    
    // stiffness
    juce::Label barrierStiffnessLabel { {}, "Barrier Stiffness (10^N)" };
    juce::Slider barrierStiffnessSlider;
    std::unique_ptr<juce::AudioProcessorValueTreeState::SliderAttachment> barrierStiffnessSliderAttachment;

    juce::Label fretStiffnessLabel { {}, "Fret Stiffness (10^N)" };
    juce::Slider fretStiffnessSlider;
    std::unique_ptr<juce::AudioProcessorValueTreeState::SliderAttachment> fretStiffnessSliderAttachment;

    juce::Label fingerStiffnessLabel { {}, "Finger Stiffness (10^N)" };
    juce::Slider fingerStiffnessSlider;
    std::unique_ptr<juce::AudioProcessorValueTreeState::SliderAttachment> fingerStiffnessSliderAttachment;

    // finger mass
    juce::Label fingerMassLabel { {}, "Finger Mass (kg)" };
    juce::Slider fingerMassSlider;
    std::unique_ptr<juce::AudioProcessorValueTreeState::SliderAttachment> fingerMassSliderAttachment;
};
