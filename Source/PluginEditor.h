/*
  ==============================================================================

    This file contains the basic framework code for a JUCE plugin editor.

  ==============================================================================
*/

#pragma once

#include <JuceHeader.h>
#include "PluginProcessor.h"
#include "DrawString.h"
#include "GuitarControlComponent.h"

//==============================================================================
/**
*/
class StringRTAudioProcessorEditor  : public juce::AudioProcessorEditor
{
public:
    StringRTAudioProcessorEditor (StringRTAudioProcessor&);
    ~StringRTAudioProcessorEditor() override;

    //==============================================================================
    void paint (juce::Graphics&) override;
    void resized() override;
    


private:
    // This reference is provided as a quick way for your editor to
    // access the processor object that created it.
    StringRTAudioProcessor& audioProcessor;
    
    GuitarControlComponent controls1;
    // GuitarControlComponent controls2;
    DrawString drawString;
    

    

    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR (StringRTAudioProcessorEditor)
};
