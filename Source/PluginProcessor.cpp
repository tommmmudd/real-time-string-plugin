/*
  ==============================================================================

    This file contains the basic framework code for a JUCE plugin processor.

  ==============================================================================
*/

#include "PluginProcessor.h"
#include "PluginEditor.h"

//==============================================================================
StringRTAudioProcessor::StringRTAudioProcessor()
#ifndef JucePlugin_PreferredChannelConfigurations
     : AudioProcessor (BusesProperties()
                     #if ! JucePlugin_IsMidiEffect
                      #if ! JucePlugin_IsSynth
                       .withInput  ("Input",  juce::AudioChannelSet::stereo(), true)
                      #endif
                       .withOutput ("Output", juce::AudioChannelSet::stereo(), true)
                     #endif
                       ),
#endif
guitarString(1), // guitarString2(2),
apvts(*this, nullptr, "ParamTree", createParameterLayout())
{
    guitarString.connectParams(apvts);
    volumeParam = apvts.getRawParameterValue("volume1");
    fretPosRangedAudioParam = apvts.getParameter("fretPosition1");
    pluckForceParam = apvts.getRawParameterValue("pluckForce1");
    pluckPosParam = apvts.getRawParameterValue("pluckPos1");
    
    volSmooth.setCurrentAndTargetValue(0);
    
}

StringRTAudioProcessor::~StringRTAudioProcessor()
{
}

//==============================================================================
const juce::String StringRTAudioProcessor::getName() const
{
    return JucePlugin_Name;
}

bool StringRTAudioProcessor::acceptsMidi() const
{
   #if JucePlugin_WantsMidiInput
    return true;
   #else
    return false;
   #endif
}

bool StringRTAudioProcessor::producesMidi() const
{
   #if JucePlugin_ProducesMidiOutput
    return true;
   #else
    return false;
   #endif
}

bool StringRTAudioProcessor::isMidiEffect() const
{
   #if JucePlugin_IsMidiEffect
    return true;
   #else
    return false;
   #endif
}

double StringRTAudioProcessor::getTailLengthSeconds() const
{
    return 0.0;
}

int StringRTAudioProcessor::getNumPrograms()
{
    return 1;   // NB: some hosts don't cope very well if you tell them there are 0 programs,
                // so this should be at least 1, even if you're not really implementing programs.
}

int StringRTAudioProcessor::getCurrentProgram()
{
    return 0;
}

void StringRTAudioProcessor::setCurrentProgram (int index)
{
}

const juce::String StringRTAudioProcessor::getProgramName (int index)
{
    return {};
}

void StringRTAudioProcessor::changeProgramName (int index, const juce::String& newName)
{
}

//==============================================================================
void StringRTAudioProcessor::prepareToPlay (double sampleRate, int samplesPerBlock)
{

    guitarString.setSampleRate(sampleRate);
    volSmooth.reset(sampleRate, 0.01);
}

void StringRTAudioProcessor::releaseResources()
{
}

#ifndef JucePlugin_PreferredChannelConfigurations
bool StringRTAudioProcessor::isBusesLayoutSupported (const BusesLayout& layouts) const
{
  #if JucePlugin_IsMidiEffect
    juce::ignoreUnused (layouts);
    return true;
  #else
    // This is the place where you check if the layout is supported.
    // In this template code we only support mono or stereo.
    // Some plugin hosts, such as certain GarageBand versions, will only
    // load plugins that support stereo bus layouts.
    if (layouts.getMainOutputChannelSet() != juce::AudioChannelSet::mono()
     && layouts.getMainOutputChannelSet() != juce::AudioChannelSet::stereo())
        return false;

    // This checks if the input layout matches the output layout
   #if ! JucePlugin_IsSynth
    if (layouts.getMainOutputChannelSet() != layouts.getMainInputChannelSet())
        return false;
   #endif

    return true;
  #endif
}
#endif

void StringRTAudioProcessor::parseMidi(juce::MidiBuffer& midi)
{
    for (const auto metadata : midi)                                                                // [9]
    {
        const auto msg = metadata.getMessage();

        if (msg.isNoteOn())
        {
            float vel = msg.getVelocity() / 128.0;
            
            // midi pitch is currently unused, but could be used to determine fret or other aspects
            float pitch = msg.getNoteNumber();
            
            float pluckPos = (*pluckPosParam);
            float pluckVel = 1.0 * (*pluckForceParam);
            pluck(1, pluckPos, vel * pluckVel, 0.001);
        }
        
        // @TODO: a particular mode where note on/offs put finger on/off?
    }

    midi.clear();
}

void StringRTAudioProcessor::processBlock (juce::AudioBuffer<double>& buffer, juce::MidiBuffer& midiMessages)
{
    juce::ScopedNoDenormals noDenormals;
    auto totalNumInputChannels  = getTotalNumInputChannels();
    auto totalNumOutputChannels = getTotalNumOutputChannels();

    for (auto i = totalNumInputChannels; i < totalNumOutputChannels; ++i)
        buffer.clear (i, 0, buffer.getNumSamples());

    // currently just using float - below
}




void StringRTAudioProcessor::processBlock (juce::AudioBuffer<float>& buffer, juce::MidiBuffer& midiMessages)
{
    juce::ScopedNoDenormals noDenormals;
    auto totalNumInputChannels  = getTotalNumInputChannels();
    auto totalNumOutputChannels = getTotalNumOutputChannels();

    for (auto i = totalNumInputChannels; i < totalNumOutputChannels; ++i)
        buffer.clear (i, 0, buffer.getNumSamples());
    
    // lazy scan for note on events
    // @TODO: make this sample accurate
    parseMidi(midiMessages);

    int numSamples = buffer.getNumSamples();
    float* left = buffer.getWritePointer(0);

    // sum all channels into single channel for audio input
    for (int i=1; i< getTotalNumInputChannels(); i++)
    {
        float* samples = buffer.getWritePointer(i);
        for (int j=0; j<numSamples; j++)
        {
            left[j] += samples[j];
        }
    }

    // add a second output channel if there is one...
    float* right;
    bool isStereo = false;
    if (getTotalNumInputChannels() > 1)
    {
        isStereo = true;
        right = buffer.getWritePointer(1);
    }

    
    guitarString.inputExternalFeedback(left, numSamples);

    // flush the left and right buffers once we've sent on the data to the strings
    buffer.clear (0, 0, numSamples);
    buffer.clear (1, 0, numSamples);
    
    // DBG( fretPosRangedAudioParam->get() );
    int fret = apvts.getParameter("fretPosition1")->convertFrom0to1 (apvts.getParameter("fretPosition1")->getValue());
    if (fret != currentFret)
    {
        currentFret = fret;
    }

    
    volSmooth.setTargetValue(*volumeParam);

    guitarString.processBlock<float>(left, right, numSamples, isStereo);

    for (int i=0; i<numSamples; i++)
    {
        float vol = volSmooth.getNextValue();
        left[i] *= vol;
        right[i] *= vol;
    }
    
    
    // could use this to automatically reinitialise
    // probably best to reset the settings that are causing the explosion in that case too though
    if (guitarString.status > 0)
        stringStatus = "Status: explosion. output set to zero. Reinitialise string.";
    else
        stringStatus = "Status: active";
    

}

//==============================================================================
bool StringRTAudioProcessor::hasEditor() const
{
    return true; // (change this to false if you choose to not supply an editor)
}

juce::AudioProcessorEditor* StringRTAudioProcessor::createEditor()
{
    return new StringRTAudioProcessorEditor (*this);
}

//==============================================================================
void StringRTAudioProcessor::getStateInformation (juce::MemoryBlock& destData)
{
    auto state = apvts.copyState();
    std::unique_ptr<juce::XmlElement> tmpXml(state.createXml());
    copyXmlToBinary(*tmpXml, destData);

}

void StringRTAudioProcessor::setStateInformation (const void* data, int sizeInBytes)
{
    std::unique_ptr<juce::XmlElement> xmlState(getXmlFromBinary(data, sizeInBytes));

    if (xmlState.get() != nullptr)
      if (xmlState->hasTagName(apvts.state.getType()))
          apvts.replaceState(juce::ValueTree::fromXml(*xmlState));
}

//==============================================================================
// This creates new instances of the plugin..
juce::AudioProcessor* JUCE_CALLTYPE createPluginFilter()
{
    return new StringRTAudioProcessor();
}
