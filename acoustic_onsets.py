# -*- coding: utf-8 -*-
"""
Article: Structural complexity of the long-term collective rhythm of 
         naturalistic conversations.
Year: 2022-2023
Author: Arodi Farrera, Caleb Rascon, Gabriel Ramos-Fernández
Source code to generate the TXT files containing acoustic onsets used as inputs
in Dynamic_complexity.R
Python version 3.8.8 (default, Apr 13 2021, 15:08:03) [MSC v.1916 64 bit (AMD64)]
Created on Tue May 18 16:44:08 2021
@author: arodi
"""


import glob
import parselmouth
import numpy as np
import os
from parselmouth.praat import call


## Parameters
timeStep_long = 0.5
timeStep_short = 0.005
minpause_long = 5.0
vowel_treshold = 0.05
int_threshold = -34

newpath = str('./results/')

def f0_array(sound, timeStep):
    
    length_samp = sound.get_number_of_samples()
    samplerate = sound.get_sampling_frequency()
    length_secs = sound.get_total_duration()
    time = np.arange(sound.xmin, sound.xmax, length_secs/length_samp)
    
    #########################
    #### ALL ABOUT PITCH ####
    #########################

    pitch = sound.to_pitch(time_step = timeStep, pitch_floor = 75, pitch_ceiling = 300)
    pitch_values = pitch.selected_array['frequency']
    min_pitch = min(pitch_values[pitch_values != 0])    
    fzero_to_data_ratio = int(length_samp/pitch_values.shape[0])
    f0_array = np.zeros(length_samp)
    
    for i in range(0, pitch_values.shape[0]):
        f0_array[i*fzero_to_data_ratio:(i+1)*fzero_to_data_ratio] = pitch_values[i]
    

    return f0_array, min_pitch, samplerate, time, length_secs, length_samp

def f0_activity_frame(f0_array, time, vowel_treshold):
    """START AND END TIME of any f0 activity"""
    vowel_starts_index = [i for i, e in enumerate(f0_array) if e > 0]
    vowel_start_frame = []
    vowel_ends_frame = []
    vowel_start_frame.append(vowel_starts_index[0])

    for (x, y) in zip(tuple(vowel_starts_index[1:]), tuple(vowel_starts_index)):
        if time[x] - time[y] > vowel_treshold:
            vowel_start_frame.append(x)
            vowel_ends_frame.append(y)
    vowel_ends_frame.append(vowel_starts_index[-1])
    
    IOIstart_vowel_timeFrame = np.array(vowel_start_frame) #when the IOI starts
    IOIend_vowel_timeFrame = np.array(vowel_ends_frame)
    
    """START AND END TIME of f0 activty"""
    vowel_x_index = [time[i] for i in IOIstart_vowel_timeFrame]
    vowel_y_index = [time[i] for i in IOIend_vowel_timeFrame-1]

    return IOIstart_vowel_timeFrame, IOIend_vowel_timeFrame, vowel_x_index, vowel_y_index

def intervention_activity_frame(f0_array, min_pitch, time, minpause):
    #################################
    #### ALL ABOUT INTERVENTIONS ####   
    #################################

    """START AND END TIME of an intervention"""

    intensity = sound.to_intensity(int(min_pitch))
    textgrid = call(intensity, "To TextGrid (silences)", int_threshold, minpause, 0.1, "silent", "sounding")
    silencetier = call(textgrid, "Extract tier", 1)
    silencetable = call(silencetier, "Down to TableOfReal", "sounding")
    npauses = call(silencetable, "Get number of rows")
    interv_start_time = []
    interv_end_time = []
 
    for i in range(1, npauses+1):
        interv_start_time.append(call(silencetable, "Get value", i, 1))
        interv_end_time.append(call(silencetable, "Get value", i, 2))

    """Only voiced (f0) intervention, either vowel or consonant (not voiceless utterances)"""
    time_round = np.round(time, 7)
    interv_starts = [round(i, 2) for i in interv_start_time]
    interv_start_frame = np.searchsorted(time_round, interv_starts)
    interv_ends = [round(i, 4) for i in interv_end_time]
    interv_end_frame = np.searchsorted(time_round, interv_ends)
    
    interv_index_bool = []
    for (x, y) in zip(tuple(interv_start_frame), tuple(interv_end_frame)):
        interv_index_bool.append(not all(f0_array[x:y] == 0))

    interv_start_frame_f0 = interv_start_frame[interv_index_bool]
    interv_end_frame_f0 = interv_end_frame[interv_index_bool]

    IOIstart_intervention_timeFrame = np.array(interv_start_frame_f0)
    IOIend_intervention_timeFrame = np.array(interv_end_frame_f0)

    """START AND END TIME of voiced interventions"""
    interv_x_index = [time[i] for i in interv_start_frame_f0]
    interv_y_index = [time[i] for i in interv_end_frame_f0-1]

    return IOIstart_intervention_timeFrame, IOIend_intervention_timeFrame, interv_x_index, interv_y_index

#for wave_file in glob.glob("*.wav", recursive = False):
#for wave_file in glob.glob("**/*.wav", recursive = True):
for wave_file in glob.glob("./mixes/EN2002b.Mix-Headset.wav", recursive = False):
    f_f0array = []
    f_f0array_short = []
    f_f0activity_start = []
    f_shortvowel = []

    
    sound = parselmouth.Sound(wave_file)
    sound = sound.resample(new_frequency = 16000)
    f_f0array = f0_array(sound, timeStep_long)
    f_interventionactivity = intervention_activity_frame(f_f0array[0], f_f0array[1], f_f0array[3], minpause_long) #f0_array, min_pitch, time

    wave_file_basename = os.path.basename(wave_file)
    print('File: ' + str(wave_file))
    print("Audio data length     = " + str(f_f0array[4]) + " secs (" + str(f_f0array[5]) + " samples)")
    print("Audio data samplerate = " + str(f_f0array[2]) + " Hz" + '\n')
   
    for (x, y) in zip(tuple(f_interventionactivity[2]), tuple(f_interventionactivity[3])):
        if (y - x) > 2.5:
            print("{:.2f}".format(y/f_f0array[4]*100)+"% :: "+str(x)+" -> "+str(y)+" : extracting")
            snd_part = sound.extract_part(from_time = x, 
                            to_time = y, preserve_times=True)
        
            f_f0array_short = f0_array(snd_part, timeStep_short)
           
            f_f0activity_start = np.array(f0_activity_frame(f_f0array_short[0], f_f0array_short[3], vowel_treshold)[2])
            vowel_onset = f_f0activity_start * f_f0array_short[2]
            f_shortvowel.append(vowel_onset)
        else:
            print("{:.2f}".format(y/f_f0array[4]*100)+"% :: "+str(x)+" -> "+str(y)+" : ---")
    print("100.0% :: Done.")
    print("---\n")
    
    
    np.savetxt(str(str(newpath) +'FRAME_shortvowel_' + str(f_f0array[5]) +'_' + str(wave_file_basename.rsplit('.', 1)[0]) + '.txt'), np.hstack(f_shortvowel), fmt='%1.0f', delimiter=",")

