# -*- coding: utf-8 -*-
"""
Article: Structural complexity of the long-term collective rhythm of 
         naturalistic conversations.
Year: 2022
Author: Arodi Farrera, Caleb Rascon, Gabriel Ramos-Fernández

Source code to generate the TXT files containing acoustic onsets used as inputs
in Dynamic_complexity.R

Created on Tue May 18 16:44:08 2021
@author: arodi
"""

import glob
import parselmouth
import numpy as np
from parselmouth.praat import call


## FOR LONG INTERVENTIONS
timeStep_long = 0.5
timeStep_short = 0.005
minpause_long =   5.0    
minpause_short = 0.3
vowel_treshold = 0.05
int_threshold = -34 

newpath = str('./')

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


def intensity_peak_activity_frame(min_pitch, samplerate):
    
    #############################
    #### ALL ABOUT INTENSITY ####   
    #############################

    intensity = sound.to_intensity(int(min_pitch))
    intensity_peaks = call(intensity, "To IntensityTier (peaks)")
    intensity_peak_table = call(intensity_peaks, "Down to TableOfReal")
    npeaks = call(intensity_peak_table, "Get number of rows")
    peak_time = []
    peak_intensity = []
    for i in range(1, npeaks+1):
        peak_time.append(call(intensity_peak_table, "Get value", i, 1)) #calls the time of the intensity peaks
        peak_intensity.append(call(intensity_peak_table, "Get value", i, 2)) #calls the intensity

    threshold = 50       
    peak_time_thres = []
    peak_intensity_thres = []
    for (x, y) in zip(tuple(peak_time), tuple(peak_intensity)):
        if y > threshold:
            peak_time_thres.append(x)
            peak_intensity_thres.append(y)

    IOI_intensity_timeFrame = np.array([int(i*samplerate) for i in peak_time_thres])
    
    return IOI_intensity_timeFrame

def get_interv_wav(sound, interv_x_index, interv_y_index, intervention_dir, wave_file):
    sound
    print("Found intervention: " + str(len(interv_x_index)) + '\n')
    for (x, y) in zip(tuple(interv_x_index), tuple(interv_y_index)):
        if x != y:
            fromTime = sound.get_time_from_frame_number(x)
            toTime = sound.get_time_from_frame_number(y)
            snd_part = sound.extract_part(from_time = fromTime, to_time = toTime)
            snd_part.save(file_path = str(intervention_dir + str(wave_file) +
                                          str(x) + '.wav'), 
                          format = 'WAV')
            print('Intervention from: ' + str(format(fromTime, '.2f')) + '(s)' + '  '+
                  'to: ' + str(format(toTime, '.2f')) + '(s)' + '\n')


for wave_file in glob.glob("*.wav", recursive = False):
    #long
    f_f0array = []
    f_interventionactivity = []
    #short
    f_f0array_short = []
    f_interventionactivity_start = []
    f_f0activity_start = []
    f_shortinterv = []
    f_shortvowel = []

    
    sound = parselmouth.Sound(wave_file)
    sound = sound.resample(new_frequency = 16000)
    f_f0array = f0_array(sound, timeStep_long)
    f_interventionactivity = intervention_activity_frame(f_f0array[0], f_f0array[1], f_f0array[3], minpause_long) #f0_array, min_pitch, time
 
    
    print('File: ' + str(wave_file))
    print("Audio data length     = " + str(f_f0array[4]) + " secs (" + str(f_f0array[5]) + " samples)")
    print("Audio data samplerate = " + str(f_f0array[2]) + " Hz" + '\n')
   
    for (x, y) in zip(tuple(f_interventionactivity[2]), tuple(f_interventionactivity[3])):
        if (y - x) > 2.5:
            snd_part = sound.extract_part(from_time = x, 
                            to_time = y, preserve_times=True)
        
            f_f0array_short = f0_array(snd_part, timeStep_short)
            f_interventionactivity_start = np.array(intervention_activity_frame(f_f0array_short[0],
                                                              f_f0array_short[1],
                                                              f_f0array_short[3],
                                                              minpause_short)[2])
            interv_onset = f_interventionactivity_start * f_f0array_short[2]
            
            
            f_f0activity_start = np.array(f0_activity_frame(f_f0array_short[0], f_f0array_short[3], vowel_treshold)[2])
            vowel_onset = f_f0activity_start * f_f0array_short[2]
            
            f_shortinterv.append(interv_onset)
            f_shortvowel.append(vowel_onset)


                  
# FRAME_longintS.txt: starts of long interventions                  
#    np.savetxt(str(str(newpath) +'FRAME_longintS_' + str(wave_file.rsplit('.', 1)[0]) + '.txt'), f_interventionactivity[0], fmt='%1.0f', delimiter=",")
# FRAME_longintE.txt: ends of long interventions 
#    np.savetxt(str(str(newpath) +'FRAME_longintE_' + str(wave_file.rsplit('.', 1)[0]) + '.txt'), f_interventionactivity[1], fmt='%1.0f', delimiter=",")
# FRAME_shortintS.txt: starts of interventions where f0 is present
#    np.savetxt(str(str(newpath) +'FRAME_shortintS_' + str(wave_file.rsplit('.', 1)[0]) + '.txt'), np.hstack(f_shortinterv), fmt='%1.0f', delimiter=",")
# FRAME_shortvowel.txt: f0 activity onsets used in dynamic_complexity.R
    np.savetxt(str(str(newpath) +'FRAME_shortvowel_' + str(f_f0array[5]) +'_' + str(wave_file.rsplit('.', 1)[0]) + '.txt'), np.hstack(f_shortvowel), fmt='%1.0f', delimiter=",")
