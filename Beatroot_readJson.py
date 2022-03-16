


#!/usr/bin/python

import numpy as np
import json
import pandas as pd

max_length_secs = 20*60
rate = 16000
minimum_intervention_length_secs = 2.5

#to read this file afterwards:
onset_file = open("EN2002b.Mix-Headset__bpm.json", "r")
audio_onset = json.load(onset_file)
onset_file.close()

max_length = int(max_length_secs*rate)
newpath = str('./')

for wave_file_basename in audio_onset:
    #print(audio_onset[wave_file_basename]["start"])
    #print(audio_onset[wave_file_basename]["end"])
    #print(audio_onset[wave_file_basename]["bpm"])
    
    data= pd.DataFrame({'start': audio_onset[wave_file_basename]["start"],
                    'end': audio_onset[wave_file_basename]["end"],
                    'bpm': audio_onset[wave_file_basename]["bpm"]})
    
    
    data.to_csv(str(str(newpath) +'BPM_' + str(wave_file_basename.rsplit('.', 1)[0]) + '.txt'),
                sep='\t')
    
