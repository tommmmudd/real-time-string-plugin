/*
  ==============================================================================

    This file contains the basic framework code for a JUCE plugin editor.

  ==============================================================================
*/

#include "PluginProcessor.h"
#include "PluginEditor.h"

//==============================================================================
StringRTAudioProcessorEditor::StringRTAudioProcessorEditor (StringRTAudioProcessor& p)
    : AudioProcessorEditor (&p), audioProcessor (p), controls1(p, p.guitarString, drawString, 1) // , controls2(p, p.guitarString2, drawString, 2)
{

    addAndMakeVisible(controls1);
    // addAndMakeVisible(controls2);
    
    addAndMakeVisible(drawString);
    
    // connect parameters from guitarString to the drawString visualisation
    drawString.setPointerAndLength (audioProcessor.guitarString.getString(), audioProcessor.guitarString.getStringLength());
    
    // for second string
    // drawString.setPointer2AndLength2 (audioProcessor.guitarString2.getString(), audioProcessor.guitarString2.getStringLength());
    
    drawString.setFretboard(audioProcessor.guitarString.getFretPositions(), audioProcessor.guitarString.getNumFrets(), audioProcessor.guitarString.getFretHeight(), audioProcessor.guitarString.getBarrierHeight());
    

    // Make sure that before the constructor has finished, you've set the
    // editor's size to whatever you need it to be.
    setSize (600, 850);
}

StringRTAudioProcessorEditor::~StringRTAudioProcessorEditor()
{
}

//==============================================================================
void StringRTAudioProcessorEditor::paint (juce::Graphics& g)
{
    // (Our component is opaque, so we must completely fill the background with a solid colour)
    //g.fillAll (getLookAndFeel().findColour (juce::ResizableWindow::backgroundColourId));
    g.fillAll (juce::Colours::darkgrey);

    g.setColour (juce::Colours::black);
    g.setFont (15.0f);
    
}

void StringRTAudioProcessorEditor::resized()
{
    int itemHeight = 30;
    controls1.setBounds(0, 0, getWidth(), getHeight());
    // controls2.setBounds(600, 0, 600, 650);
    
    auto stringArea = getLocalBounds().removeFromBottom(itemHeight*7).reduced(10);
    
    drawString.setBounds(stringArea.removeFromTop(itemHeight*7));
}
