/*
  ==============================================================================

    GuitarControlComponent.cpp
    Created: 11 Feb 2025 3:33:40pm
    Author:  Tom Mudd

  ==============================================================================
*/

#include "GuitarControlComponent.h"


GuitarControlComponent::GuitarControlComponent(StringRTAudioProcessor& p, GuitarString& gs, DrawString& ds, int ind) : audioProcessor(p), guitarString(gs), drawString(ds), index(ind)
{
    
    i_s = std::to_string(index);
    
    addAndMakeVisible(resetButton);
    resetButton.onClick = [this] { audioProcessor.resetString(index); };   // [8]
    resetButton.setButtonText ("Reinitialise string");
    
    addAndMakeVisible(pluckButton);
    pluckButton.onClick = [this] { audioProcessor.pluck(index, *pluckPosParam, *pluckForceParam, 0.001); };   // [8]
    pluckButton.setButtonText ("Pluck");
    
    
    // Real Time Params
    
    addAndMakeVisible(fingerPosLabel);
    addAndMakeVisible(fingerPosSlider);
    fingerPosSlider.setRange(0.0f, 0.85f, 0.001f);
    //fingerPosSlider.setTextBoxStyle (juce::Slider::NoTextBox, false, 60, fingerPosSlider.getTextBoxHeight());
    fingerPosSlider.setValue(0.25f);
    fingerPosSliderAttachment = std::make_unique<juce::AudioProcessorValueTreeState::SliderAttachment>( audioProcessor.apvts, "fingerPos"+i_s, fingerPosSlider);
    
    

    
    addAndMakeVisible(fretPositionLabel);
    addAndMakeVisible(fretPositionSlider);
    fretPositionSlider.setRange(1, 20, 1);
    fretPositionSlider.setValue(1);
    fretPositionSliderAttachment = std::make_unique<juce::AudioProcessorValueTreeState::SliderAttachment>( audioProcessor.apvts, "fretPosition"+i_s, fretPositionSlider);
    fretPositionSlider.onValueChange = [this] { putFingerAtFret(fretPositionSlider.getValue()); };

    
    addAndMakeVisible(fingerForceLabel);
    addAndMakeVisible(fingerForceSlider);
    fingerForceSlider.setRange(-0.05f, 1.0f, 0.0001f);     // optional 3rd argument is interval
    //fingerForceSlider.setTextBoxStyle (juce::Slider::NoTextBox, false, 60, fingerForceSlider.getTextBoxHeight());
    fingerForceSlider.setValue(0.0f);
    fingerForceSliderAttachment = std::make_unique<juce::AudioProcessorValueTreeState::SliderAttachment>( audioProcessor.apvts, "fingerForce"+i_s, fingerForceSlider);
    
    addAndMakeVisible(outputPos1Label);
    addAndMakeVisible(outputPos1Slider);
    outputPos1Slider.setRange(0.001f, 0.999f, 0.001f);
    outputPos1Slider.setValue(0.72f);
    outputPos1SliderAttachment = std::make_unique<juce::AudioProcessorValueTreeState::SliderAttachment>( audioProcessor.apvts, "outputPos1"+i_s, outputPos1Slider);
    
    addAndMakeVisible(outputPos2Label);
    addAndMakeVisible(outputPos2Slider);
    outputPos2Slider.setRange(0.001f, 0.999f, 0.001f);
    outputPos2Slider.setValue(0.72f);
    outputPos2SliderAttachment = std::make_unique<juce::AudioProcessorValueTreeState::SliderAttachment>( audioProcessor.apvts, "outputPos2"+i_s, outputPos2Slider);

    
    addAndMakeVisible(externalFeedbackGainLabel);
    addAndMakeVisible(externalFeedbackGainSlider);
    // externalFeedbackGainSlider.setRange(0.01f, 5.0f, 0.01f);
    // externalFeedbackGainSlider.setValue(0.2f);
    externalFeedbackSliderAttachment = std::make_unique<juce::AudioProcessorValueTreeState::SliderAttachment>( audioProcessor.apvts, "externalFeedbackGain"+i_s, externalFeedbackGainSlider);
    
    addAndMakeVisible(internalFeedbackGainLabel);
    addAndMakeVisible(internalFeedbackGainSlider);
    // internalFeedbackGainSlider.setRange(0.01f, 5.0f, 0.01f);
    // internalFeedbackGainSlider.setValue(0.2f);
    internalFeedbackGainSliderAttachment = std::make_unique<juce::AudioProcessorValueTreeState::SliderAttachment>( audioProcessor.apvts, "feedbackGain"+i_s, internalFeedbackGainSlider);
    
    addAndMakeVisible(feedbackFilterLabel);
    addAndMakeVisible(feedbackFilterSlider);
    feedbackFilterSlider.setRange(0.01f, 5.0f, 0.01f);
    feedbackFilterSlider.setValue(0.2f);
    feedbackFilterSliderAttachment = std::make_unique<juce::AudioProcessorValueTreeState::SliderAttachment>( audioProcessor.apvts, "feedbackFilter"+i_s, feedbackFilterSlider);
    
    addAndMakeVisible(feedbackPosLabel);
    addAndMakeVisible(feedbackPosSlider);
    feedbackPosSlider.setRange(0.01f, 5.0f, 0.01f);
    feedbackPosSlider.setValue(0.2f);
    feedbackPosSliderAttachment = std::make_unique<juce::AudioProcessorValueTreeState::SliderAttachment>( audioProcessor.apvts, "feedbackPos"+i_s, feedbackPosSlider);
    
    
    
    
    addAndMakeVisible(pluckForceLabel);
    addAndMakeVisible(pluckForceSlider);
    pluckForceSlider.setRange(0.01f, 5.0f, 0.01f);
    pluckForceSlider.setValue(0.2f);
    pluckForceSliderAttachment = std::make_unique<juce::AudioProcessorValueTreeState::SliderAttachment>( audioProcessor.apvts, "pluckForce"+i_s, pluckForceSlider);
    
    
    addAndMakeVisible(pluckPosLabel);
    addAndMakeVisible(pluckPosSlider);
    pluckPosSlider.setRange(0.01f, 1.0f, 0.01f);
    pluckPosSlider.setValue(0.88f);
    pluckPosSliderAttachment = std::make_unique<juce::AudioProcessorValueTreeState::SliderAttachment>( audioProcessor.apvts, "pluckPos"+i_s, pluckPosSlider);
    
    
    
    
    addAndMakeVisible(kirchButton);
    kirchButton.onClick = [this] { setFlag (&kirchButton, "kirch"); };
    kirchButton.setToggleState(true, juce::NotificationType::dontSendNotification);
    addAndMakeVisible(fretsButton);
    fretsButton.onClick = [this] { setFlag (&fretsButton, "frets"); };
    fretsButton.setToggleState(true, juce::NotificationType::dontSendNotification);
    addAndMakeVisible(barrierButton);
    barrierButton.onClick = [this] { setFlag (&barrierButton, "barrier"); };
    barrierButton.setToggleState(true, juce::NotificationType::dontSendNotification);
    // kirchButton.setClickingTogglesState (true);
    addAndMakeVisible(fingerButton);
    fingerButton.onClick = [this] { setFlag (&fingerButton, "finger"); };
    fingerButton.setToggleState(true, juce::NotificationType::dontSendNotification);
    // kirchButton.setClickingTogglesState (true);
    
    kcAttachment = std::make_unique<juce::AudioProcessorValueTreeState::ButtonAttachment>( audioProcessor.apvts, "kcToggle"+i_s, kirchButton);
    barrierAttachment = std::make_unique<juce::AudioProcessorValueTreeState::ButtonAttachment>( audioProcessor.apvts, "barrierToggle"+i_s, barrierButton);
    fretsAttachment = std::make_unique<juce::AudioProcessorValueTreeState::ButtonAttachment>( audioProcessor.apvts, "fretsToggle"+i_s, fretsButton);
    fingerAttachment = std::make_unique<juce::AudioProcessorValueTreeState::ButtonAttachment>( audioProcessor.apvts, "fingerToggle"+i_s, fingerButton);
    
    
    addAndMakeVisible(volumeLabel);
    addAndMakeVisible(volumeSlider);
    volumeSlider.setRange(0.0f, 1.0f, 0.0001f);
    volumeSlider.setValue(0.5f);
    volumeSliderAttachment = std::make_unique<juce::AudioProcessorValueTreeState::SliderAttachment>( audioProcessor.apvts, "volume"+i_s, volumeSlider);
    
    // ok to do this here? multiple pointers?
    pluckForceParam = audioProcessor.apvts.getRawParameterValue("pluckForce"+i_s);
    pluckPosParam   = audioProcessor.apvts.getRawParameterValue("pluckPos"+i_s);
    
    addAndMakeVisible(statusLabel);
    
    
    
    // STRING PARAMS
    addAndMakeVisible(stringLengthLabel);
    addAndMakeVisible(stringLengthSlider);
    stringLengthSlider.setRange(0.3f, 2.4f, 0.01f);
    stringLengthSlider.setValue(0.68f);
    stringLengthSliderAttachment = std::make_unique<juce::AudioProcessorValueTreeState::SliderAttachment>( audioProcessor.apvts, "stringLength"+i_s, stringLengthSlider);

    addAndMakeVisible(stringTensionLabel);
    addAndMakeVisible(stringTensionSlider);
    stringTensionSlider.setRange(4.0f, 60.0f, 0.1f);
    stringTensionSlider.setValue(12.1f);
    stringTensionSliderAttachment = std::make_unique<juce::AudioProcessorValueTreeState::SliderAttachment>( audioProcessor.apvts, "stringTension"+i_s, stringTensionSlider);
    
    addAndMakeVisible(stringRadiusLabel);
    addAndMakeVisible(stringRadiusSlider);
    stringRadiusSlider.setRange(5e-05f, 0.001f, 0.00001f);
    stringRadiusSlider.setValue(0.0002f);
    stringRadiusSliderAttachment = std::make_unique<juce::AudioProcessorValueTreeState::SliderAttachment>( audioProcessor.apvts, "stringRadius"+i_s, stringRadiusSlider);

    addAndMakeVisible(stringT60_0Label);
    addAndMakeVisible(stringT60_0Slider);
    stringT60_0Slider.setRange(2.0f, 50.0f, 0.1f);
    stringT60_0Slider.setValue(15.0f);
    stringT60_0SliderAttachment = std::make_unique<juce::AudioProcessorValueTreeState::SliderAttachment>( audioProcessor.apvts, "stringT60_0"+i_s, stringT60_0Slider);

    addAndMakeVisible(stringT60_1Label);
    addAndMakeVisible(stringT60_1Slider);
    stringT60_1Slider.setRange(1.0f, 48.0f, 0.1f);
    stringT60_1Slider.setValue(5.0f);
    stringT60_1SliderAttachment = std::make_unique<juce::AudioProcessorValueTreeState::SliderAttachment>( audioProcessor.apvts, "stringT60_1"+i_s, stringT60_1Slider);


    
    // Fret height
    addAndMakeVisible(fretHeightLabel);
    addAndMakeVisible(fretHeightSlider);
    fretHeightSlider.setRange(-0.01f, 0.0f, 0.00005f);
    fretHeightSlider.setValue(-0.0005f);
    fretHeightSliderAttachment = std::make_unique<juce::AudioProcessorValueTreeState::SliderAttachment>( audioProcessor.apvts, "fretHeight"+i_s, fretHeightSlider);
    
    
    // barrier height
    addAndMakeVisible(barrierHeightLabel);
    addAndMakeVisible(barrierHeightSlider);
    barrierHeightSlider.setRange(-0.02f, 0.0f, 0.0001f);
    barrierHeightSlider.setValue(-0.001f);
    barrierHeightSliderAttachment = std::make_unique<juce::AudioProcessorValueTreeState::SliderAttachment>( audioProcessor.apvts, "barrierHeight"+i_s, barrierHeightSlider);
    
    
    
    
    

    // ALPHAS
    addAndMakeVisible(alphaBarrierLabel);
    addAndMakeVisible(alphaBarrierSlider);
    alphaBarrierSlider.setRange(1.0, 3.5f, 0.0001f);
    alphaBarrierSlider.setValue(2.3f);
    alphaBarrierSliderAttachment = std::make_unique<juce::AudioProcessorValueTreeState::SliderAttachment>( audioProcessor.apvts, "barrierAlpha"+i_s, alphaBarrierSlider);
    
    addAndMakeVisible(alphaFretLabel);
    addAndMakeVisible(alphaFretSlider);
    alphaFretSlider.setRange(1.0, 3.5f, 0.001f);
    alphaFretSlider.setValue(2.3f);
    alphaFretSliderAttachment = std::make_unique<juce::AudioProcessorValueTreeState::SliderAttachment>( audioProcessor.apvts, "fretAlpha"+i_s, alphaFretSlider);
    
    addAndMakeVisible(alphaFingerLabel);
    addAndMakeVisible(alphaFingerSlider);
    alphaFingerSlider.setRange(1.0, 3.5f, 0.001f);
    alphaFingerSlider.setValue(1.3f);
    alphaFingerSliderAttachment = std::make_unique<juce::AudioProcessorValueTreeState::SliderAttachment>( audioProcessor.apvts, "fingerAlpha"+i_s, alphaFingerSlider);
    
    
    
    // stiffness
    addAndMakeVisible(barrierStiffnessLabel);
    addAndMakeVisible(barrierStiffnessSlider);
    barrierStiffnessSlider.setRange(2, 20, 1);
    barrierStiffnessSlider.setValue(15);
    barrierStiffnessSliderAttachment = std::make_unique<juce::AudioProcessorValueTreeState::SliderAttachment>( audioProcessor.apvts, "barrierStiffness"+i_s, barrierStiffnessSlider);

    addAndMakeVisible(fretStiffnessLabel);
    addAndMakeVisible(fretStiffnessSlider);
    fretStiffnessSlider.setRange(2, 20, 1);
    fretStiffnessSlider.setValue(15);
    fretStiffnessSliderAttachment = std::make_unique<juce::AudioProcessorValueTreeState::SliderAttachment>( audioProcessor.apvts, "fretStiffness"+i_s, fretStiffnessSlider);

    addAndMakeVisible(fingerStiffnessLabel);
    addAndMakeVisible(fingerStiffnessSlider);
    fingerStiffnessSlider.setRange(2, 20, 1);
    fingerStiffnessSlider.setValue(10);
    fingerStiffnessSliderAttachment = std::make_unique<juce::AudioProcessorValueTreeState::SliderAttachment>( audioProcessor.apvts, "fingerStiffness"+i_s, fingerStiffnessSlider);

    

    // finger mass
    addAndMakeVisible(fingerMassLabel);
    addAndMakeVisible(fingerMassSlider);
    fingerMassSlider.setRange(0.0001f, 4.0f, 0.0001f);
    fingerMassSlider.setValue(0.005f);
    fingerMassSliderAttachment = std::make_unique<juce::AudioProcessorValueTreeState::SliderAttachment>( audioProcessor.apvts, "fingerMass"+i_s, fingerMassSlider);
}


void GuitarControlComponent::paint(juce::Graphics& g)
{
    g.fillAll (juce::Colours::darkgrey);

    g.setColour (juce::Colours::black);
    g.setFont (15.0f);
    
    if (index == 1)
    {
        statusLabel.setText(audioProcessor.stringStatus, juce::NotificationType::dontSendNotification);
    }
    else
    {
        // statusLabel.setText(audioProcessor.stringStatus2, juce::NotificationType::dontSendNotification);
    }
}


void GuitarControlComponent::resized()
{
    // left side (real time)
    auto leftArea = getLocalBounds().removeFromLeft(getWidth()/2).reduced(10);
    int itemHeight = 22;
    
    resetButton.setBounds(leftArea.removeFromTop(itemHeight));
    pluckButton.setBounds(leftArea.removeFromTop(itemHeight));
    
    fingerPosLabel.setBounds(leftArea.removeFromTop(itemHeight));
    fingerPosSlider.setBounds(leftArea.removeFromTop(itemHeight));
    
    fretPositionLabel.setBounds(leftArea.removeFromTop(itemHeight));
    fretPositionSlider.setBounds(leftArea.removeFromTop(itemHeight));
    
    fingerForceLabel.setBounds(leftArea.removeFromTop(itemHeight));
    fingerForceSlider.setBounds(leftArea.removeFromTop(itemHeight));
    outputPos1Label.setBounds(leftArea.removeFromTop(itemHeight));
    outputPos1Slider.setBounds(leftArea.removeFromTop(itemHeight));
    outputPos2Label.setBounds(leftArea.removeFromTop(itemHeight));
    outputPos2Slider.setBounds(leftArea.removeFromTop(itemHeight));
    
    externalFeedbackGainLabel.setBounds(leftArea.removeFromTop(itemHeight));
    externalFeedbackGainSlider.setBounds(leftArea.removeFromTop(itemHeight));
    
    internalFeedbackGainLabel.setBounds(leftArea.removeFromTop(itemHeight));
    internalFeedbackGainSlider.setBounds(leftArea.removeFromTop(itemHeight));
    feedbackPosLabel.setBounds(leftArea.removeFromTop(itemHeight));
    feedbackPosSlider.setBounds(leftArea.removeFromTop(itemHeight));
    feedbackFilterLabel.setBounds(leftArea.removeFromTop(itemHeight));
    feedbackFilterSlider.setBounds(leftArea.removeFromTop(itemHeight));
    
    pluckForceLabel.setBounds(leftArea.removeFromTop(itemHeight));
    pluckForceSlider.setBounds(leftArea.removeFromTop(itemHeight));
    pluckPosLabel.setBounds(leftArea.removeFromTop(itemHeight));
    pluckPosSlider.setBounds(leftArea.removeFromTop(itemHeight));
    
    // buttons as 2x2 grid
    auto buttonArea = leftArea.removeFromTop(itemHeight * 1.5);
    kirchButton.setBounds(buttonArea);
    fretsButton.setBounds(buttonArea.getX() + 100, buttonArea.getY(), buttonArea.getWidth(), buttonArea.getHeight());
    
    buttonArea = leftArea.removeFromTop(itemHeight);
    barrierButton.setBounds(buttonArea);
    fingerButton.setBounds(buttonArea.getX() + 100, buttonArea.getY(), buttonArea.getWidth(), buttonArea.getHeight());
    
    
    
    itemHeight = 20;
    auto rightArea = getLocalBounds().removeFromRight(300).reduced(10);
    
    stringLengthLabel.setBounds(rightArea.removeFromTop(itemHeight));
    stringLengthSlider.setBounds(rightArea.removeFromTop(itemHeight));
    stringTensionLabel.setBounds(rightArea.removeFromTop(itemHeight));
    stringTensionSlider.setBounds(rightArea.removeFromTop(itemHeight));
    stringRadiusLabel.setBounds(rightArea.removeFromTop(itemHeight));
    stringRadiusSlider.setBounds(rightArea.removeFromTop(itemHeight));
    stringT60_0Label.setBounds(rightArea.removeFromTop(itemHeight));
    stringT60_0Slider.setBounds(rightArea.removeFromTop(itemHeight));
    stringT60_1Label.setBounds(rightArea.removeFromTop(itemHeight));
    stringT60_1Slider.setBounds(rightArea.removeFromTop(itemHeight));
    
    fretHeightLabel.setBounds(rightArea.removeFromTop(itemHeight));
    fretHeightSlider.setBounds(rightArea.removeFromTop(itemHeight));
    barrierHeightLabel.setBounds(rightArea.removeFromTop(itemHeight));
    barrierHeightSlider.setBounds(rightArea.removeFromTop(itemHeight));
    
    alphaBarrierLabel.setBounds(rightArea.removeFromTop(itemHeight));
    alphaBarrierSlider.setBounds(rightArea.removeFromTop(itemHeight));
    alphaFretLabel.setBounds(rightArea.removeFromTop(itemHeight));
    alphaFretSlider.setBounds(rightArea.removeFromTop(itemHeight));
    alphaFingerLabel.setBounds(rightArea.removeFromTop(itemHeight));
    alphaFingerSlider.setBounds(rightArea.removeFromTop(itemHeight));
    
    barrierStiffnessLabel.setBounds(rightArea.removeFromTop(itemHeight));
    barrierStiffnessSlider.setBounds(rightArea.removeFromTop(itemHeight));
    fretStiffnessLabel.setBounds(rightArea.removeFromTop(itemHeight));
    fretStiffnessSlider.setBounds(rightArea.removeFromTop(itemHeight));
    fingerStiffnessLabel.setBounds(rightArea.removeFromTop(itemHeight));
    fingerStiffnessSlider.setBounds(rightArea.removeFromTop(itemHeight));
    
    fingerMassLabel.setBounds(rightArea.removeFromTop(itemHeight));
    fingerMassSlider.setBounds(rightArea.removeFromTop(itemHeight));
    
    statusLabel.setBounds(rightArea.removeFromTop(itemHeight));
    
    
    volumeLabel.setBounds(leftArea.removeFromTop(itemHeight*1.5));
    auto volSliderArea = leftArea.removeFromTop(itemHeight);
    volumeSlider.setBounds(volSliderArea.getX(), volSliderArea.getY(), getWidth() - 20, volSliderArea.getHeight());
}
