/*
  ==============================================================================

    This file contains the basic framework code for a JUCE plugin processor.

  ==============================================================================
*/

#pragma once

#include <JuceHeader.h>
#include "GuitarString.h"

//==============================================================================
/**
*/
class StringRTAudioProcessor  : public juce::AudioProcessor
{
public:
    //==============================================================================
    StringRTAudioProcessor();
    ~StringRTAudioProcessor() override;

    //==============================================================================
    void prepareToPlay (double sampleRate, int samplesPerBlock) override;
    void releaseResources() override;

   #ifndef JucePlugin_PreferredChannelConfigurations
    bool isBusesLayoutSupported (const BusesLayout& layouts) const override;
   #endif

    void processBlock (juce::AudioBuffer<double>&, juce::MidiBuffer&) override;
    
    void processBlock (juce::AudioBuffer<float>&, juce::MidiBuffer&) override;

    //==============================================================================
    juce::AudioProcessorEditor* createEditor() override;
    bool hasEditor() const override;

    //==============================================================================
    const juce::String getName() const override;

    bool acceptsMidi() const override;
    bool producesMidi() const override;
    bool isMidiEffect() const override;
    double getTailLengthSeconds() const override;

    //==============================================================================
    int getNumPrograms() override;
    int getCurrentProgram() override;
    void setCurrentProgram (int index) override;
    const juce::String getProgramName (int index) override;
    void changeProgramName (int index, const juce::String& newName) override;

    //==============================================================================
    void getStateInformation (juce::MemoryBlock& destData) override;
    void setStateInformation (const void* data, int sizeInBytes) override;
    
    /// pluck a particular string at a particular point (currently only one string, so whichString unused)
    void pluck(int whichString, float pos, float force, float dur)
    {
        guitarString.pluck(pos, force, dur);
    }
    
    /// if it goes awry, can reset it to go again
    void resetString(int whichString)
    {
        guitarString.setup();

    }
    
    juce::AudioProcessorValueTreeState apvts;
    
    // the actual string DSP
    GuitarString guitarString;
    
    juce::RangedAudioParameter* fretPosRangedAudioParam;
    
    juce::String stringStatus {"Status: active"};

private:

    
    
    std::atomic<float>* volumeParam;
    std::atomic<float>* pluckForceParam;
    std::atomic<float>* pluckPosParam;
    
    juce::SmoothedValue<float> volSmooth;
    
    int currentFret = 1;
    
    void parseMidi(juce::MidiBuffer& midiMessages);
    
    
    /// initialise plugin parameters
    juce::AudioProcessorValueTreeState::ParameterLayout createParameterLayout()
    {
        juce::AudioProcessorValueTreeState::ParameterLayout layout;
        
        juce::NormalisableRange<float> normRangePos (0.0f, 0.85f, 0.0001f, 1);
        juce::NormalisableRange<float> normRangeForce (0.0f, 1.0f, 0.0001f, 1);
        
        
        // Loop for multiple strings (currently just one)
        for (int i=1; i<2; i++)
        {
            std::string i_s = std::to_string(i);
            
            // REAL TIME PARAMS
            layout.add(std::make_unique<juce::AudioParameterFloat>(juce::ParameterID("fingerPos"+i_s, 1), "Finger Pos", normRangePos, 0.25f));
            
            layout.add(std::make_unique<juce::AudioParameterInt>(juce::ParameterID("fretPosition"+i_s, 1), "Finger at Fret", 1, 20, 1));
            
            
            layout.add(std::make_unique<juce::AudioParameterFloat>(juce::ParameterID("fingerForce"+i_s, 1), "Finger Force", normRangeForce, 0.0f));
            
            layout.add(std::make_unique<juce::AudioParameterFloat>(juce::ParameterID("outputPos1"+i_s, 1), "Output Pos 1", juce::NormalisableRange<float>(0.001f, 0.999f, 0.0001f, 1), 0.72f));
            layout.add(std::make_unique<juce::AudioParameterFloat>(juce::ParameterID("outputPos2"+i_s, 1), "Output Pos 2", juce::NormalisableRange<float>(0.001f, 0.999f, 0.0001f, 1), 0.89f));
            
            layout.add(std::make_unique<juce::AudioParameterFloat>(juce::ParameterID("externalFeedbackGain"+i_s, 1), "External Feedback Gain", juce::NormalisableRange<float>(0.0f, 12.0f, 0.01f, 0.5), 0.0f));
            
            layout.add(std::make_unique<juce::AudioParameterFloat>(juce::ParameterID("feedbackGain"+i_s, 1), "Feedback Gain", juce::NormalisableRange<float>(0.0f, 12.0f, 0.01f, 0.5), 0.0f));
            
            layout.add(std::make_unique<juce::AudioParameterFloat>(juce::ParameterID("feedbackPos"+i_s, 1), "Feedback Pos", juce::NormalisableRange<float>(0.001f, 0.999f, 0.0001f, 1), 0.5f));
            
            layout.add(std::make_unique<juce::AudioParameterFloat>(juce::ParameterID("feedbackFilter"+i_s, 1), "Feedback Filter", juce::NormalisableRange<float>(100.0f, 15000.0f, 1.0f, 0.5), 3500.0f));
            
            layout.add(std::make_unique<juce::AudioParameterFloat>(juce::ParameterID("pluckForce"+i_s, 1), "Pluck Force", juce::NormalisableRange<float>(0.01f, 5.0f, 0.01f, 0.5), 0.2f));
            
            layout.add(std::make_unique<juce::AudioParameterFloat>(juce::ParameterID("pluckPos"+i_s, 1), "Pluck Position", juce::NormalisableRange<float>(0.001f, 0.999f, 0.0001f, 1), 0.88f));
            
            layout.add(std::make_unique<juce::AudioParameterFloat>(juce::ParameterID("volume"+i_s, 1), "Volume", juce::NormalisableRange<float>(0.0f, 1.0f, 0.01f, 1.0), 0.5f));
            
            
            
            // STRING PARAMS
            layout.add(std::make_unique<juce::AudioParameterFloat>(juce::ParameterID("stringLength"+i_s, 1), "String Length", 0.3f, 4.8f, 0.68f));
            layout.add(std::make_unique<juce::AudioParameterFloat>(juce::ParameterID("stringTension"+i_s, 1), "String Tension", juce::NormalisableRange<float>(4.0f, 60.0f, 0.01f, 1), 12.1f));
            layout.add(std::make_unique<juce::AudioParameterFloat>(juce::ParameterID("stringRadius"+i_s, 1), "String Radius", juce::NormalisableRange<float>(5e-05f, 0.001f, 0.00001f, 1), 0.0002f));
            layout.add(std::make_unique<juce::AudioParameterFloat>(juce::ParameterID("stringT60_0"+i_s, 1), "String Low Decay", juce::NormalisableRange<float>(2.0f, 50.0f, 0.1f, 1), 15.0f));
            layout.add(std::make_unique<juce::AudioParameterFloat>(juce::ParameterID("stringT60_1"+i_s, 1), "String High Decay", juce::NormalisableRange<float>(1.0f, 48.0f, 0.1f, 1), 5.0f));
            
            
            // FRET PARAMS
            layout.add(std::make_unique<juce::AudioParameterFloat>(juce::ParameterID("fretHeight"+i_s, 1), "Fret height", juce::NormalisableRange<float>(-0.01f, 0.0f, 0.0001f, 2), -0.0005f));
            
            layout.add(std::make_unique<juce::AudioParameterFloat>(juce::ParameterID("barrierHeight"+i_s, 1), "Barrier height", juce::NormalisableRange<float>(-0.02f, 0.0f, 0.0001f, 2), -0.001f));
            
            
            layout.add(std::make_unique<juce::AudioParameterFloat>(juce::ParameterID("barrierAlpha"+i_s, 1), "Barrier Alpha", 1.0f, 3.5f, 2.3f));
            layout.add(std::make_unique<juce::AudioParameterFloat>(juce::ParameterID("fretAlpha"+i_s, 1), "Fret Alpha", 1.0f, 3.5f, 2.3f));
            layout.add(std::make_unique<juce::AudioParameterFloat>(juce::ParameterID("fingerAlpha"+i_s, 1), "Finger Alpha", 1.0f, 3.5f, 1.3f));
            
            
            layout.add(std::make_unique<juce::AudioParameterFloat>(juce::ParameterID("barrierStiffness"+i_s, 1), "Barrier Stiffness (10^N)", juce::NormalisableRange<float>(2, 20, 1, 1), 15));
            
            layout.add(std::make_unique<juce::AudioParameterFloat>(juce::ParameterID("fretStiffness"+i_s, 1), "Fret Stiffness (10^N)", juce::NormalisableRange<float>(2, 20, 1, 1), 15));
            
            layout.add(std::make_unique<juce::AudioParameterFloat>(juce::ParameterID("fingerStiffness"+i_s, 1), "Finger Stiffness (10^N)", juce::NormalisableRange<float>(2, 20, 1, 1), 10));
            
            layout.add(std::make_unique<juce::AudioParameterFloat>(juce::ParameterID("fingerMass"+i_s, 1), "Finger Mass (kg)", juce::NormalisableRange<float>(0.0001f, 4.0f, 0.0001f, 0.5), 0.005f));
            
            
            layout.add(std::make_unique<juce::AudioParameterBool>(juce::ParameterID("kcToggle"+i_s, 1), "Kirchoff", true));
            layout.add(std::make_unique<juce::AudioParameterBool>(juce::ParameterID("barrierToggle"+i_s, 1), "Barrier", true));
            layout.add(std::make_unique<juce::AudioParameterBool>(juce::ParameterID("fretsToggle"+i_s, 1), "Frets", true));
            layout.add(std::make_unique<juce::AudioParameterBool>(juce::ParameterID("fingerToggle"+i_s, 1), "Finger", true));
            
        }
        
        
        
        return layout;
    }
    
    //==============================================================================
    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR (StringRTAudioProcessor)
};
