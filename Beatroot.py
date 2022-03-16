"""
Article: Structural complexity of the long-term collective rhythm of 
         naturalistic conversations.
Year: 2022
Author: Arodi Farrera, Caleb Rascon, Gabriel Ramos-Fernández

Source code to generate the TXT files containing BPM values used as inputs
in Dynamic_complexity.R

@author: caleb
"""


#!/usr/bin/python3

# full paper: https://www.tandfonline.com/doi/abs/10.1076/jnmr.30.1.39.7119
# implementation 'borrowed' from:
#    https://github.com/JZacaroli/Beat-Tracking/blob/master/BeatDetection.ipynb
#
# requires Parselmouth:
#    in Ubuntu 20.04: pip3 install praat-parselmouth
#    in Ubuntu 18.04, run with python3 and:
#                     pip3 install cmake
#                     pip3 install scikit-build
#                     pip3 install praat-parselmouth
#                     pip3 install soundfile

import os
import sys
import glob
import parselmouth
import numpy as np
import json
from parselmouth.praat import call
import soundfile

def f0_array(sound, timeStep = 0.005):
    
    length_samp = sound.get_number_of_samples()
    samplerate = sound.get_sampling_frequency()
    length_secs = sound.get_total_duration()
    time = np.arange(sound.xmin, sound.xmax, length_secs/length_samp)
    
    #########################
    #### ALL ABOUT PITCH ####
    #########################

# Obtiene valores de f0 (sound.to_pitch, pitch.selected_array)
    pitch = sound.to_pitch(time_step = timeStep, pitch_floor = 60, pitch_ceiling = 600)
    pitch_values = pitch.selected_array['frequency']
    min_pitch = min(pitch_values[pitch_values != 0])

# La clase de los objetos que regresan esas funciones son muy especificas
# asi que los pasos siguientes son para obtener un vector (f0_array) 
# con longitud = numero de samples y donde los valores de f0 esten en los frames correspondientes 
    
    fzero_to_data_ratio = int(length_samp/pitch_values.shape[0])
    f0_array = np.zeros(length_samp)
    
    for i in range(0, pitch_values.shape[0]):
        f0_array[i*fzero_to_data_ratio:(i+1)*fzero_to_data_ratio] = pitch_values[i]
    
# Regresa f0_array y otros datos que se ocupan en otras funciones
# length_sec = duracion en segundos; length_samp = numero de muestras
    return f0_array, min_pitch, samplerate, time, length_secs, length_samp



    
def f0_activity_frame(f0_array, time, vowel_treshold = 0.05):
    """START AND END TIME of any f0 activity"""
# Obtiene los frames en donde el f0 es mayor a 0 (ahora que lo pienso parece redundante
# pero estoy segura que tiene una explicacion logica)
    vowel_starts_index = [i for i, e in enumerate(f0_array) if e > 0]
    vowel_start_frame = []
    vowel_ends_frame = []
    vowel_start_frame.append(vowel_starts_index[0])

# si el tiempo que pasa entre la actividad Y y X supera un umbral (vowel_threshold)
# seguro es una nueva vocal/silaba
# se guarda el frame de esa nueva vocal

    for (x, y) in zip(tuple(vowel_starts_index[1:]), tuple(vowel_starts_index)):
        if time[x] - time[y] > vowel_treshold:
            vowel_start_frame.append(x)
            vowel_ends_frame.append(y)
    vowel_ends_frame.append(vowel_starts_index[-1])
    
    IOIstart_vowel_timeFrame = np.array(vowel_start_frame) #when the IOI starts
    IOIend_vowel_timeFrame = np.array(vowel_ends_frame)

# Regresa vectores con los valores de inicio (IOIstart) y final (IOIend)
# de una vocal (i.e. actividad de f0 delimitada por espacios sin actividad de f0 
# > vowel_treshold  -me falto una h en threshold)
    return IOIstart_vowel_timeFrame, IOIend_vowel_timeFrame


def intervention_activity_frame(sound, f0_array, min_pitch, time, int_threshold =  -34, minpause = 0.3):
    #################################
    #### ALL ABOUT INTERVENTIONS ####   
    #################################

    """START AND END TIME of an intervention"""
# Funciones para obtener actividad acustica de Praat
# (to_intensity): buscar sound:to intensity
# (intensity): buscar Intensity: to textgrid(silences)

    intensity = sound.to_intensity(int(min_pitch))
    textgrid = call(intensity, "To TextGrid (silences)", int_threshold, minpause, 0.1, "silent", "sounding")
    silencetier = call(textgrid, "Extract tier", 1)
# Pasos solo para volver el objeto a una tabla de datos numericos
    silencetable = call(silencetier, "Down to TableOfReal", "sounding")
    npauses = call(silencetable, "Get number of rows")
    interv_start_time = []
    interv_end_time = []
 
# tabla con el inicio y final (s) de las intervenciones    
    for i in range(1, npauses+1):
        interv_start_time.append(call(silencetable, "Get value", i, 1))
        interv_end_time.append(call(silencetable, "Get value", i, 2))

# Pasos para saber si en el rango temporal de las intervenciones 
# existen valores de f0

# obtiene el frame equivalente a los datos temporales obtenidos en el paso anterior
    """Only voiced (f0) intervention, either vowel or consonant (not voiceless utterances)"""
    time_round = np.round(time, 7)
    interv_starts = [round(i, 2) for i in interv_start_time]
    interv_start_frame = np.searchsorted(time_round, interv_starts)
    interv_ends = [round(i, 4) for i in interv_end_time]
    interv_end_frame = np.searchsorted(time_round, interv_ends)
    
# vector logico sobre si existen valores de f0
    interv_index_bool = []
    for (x, y) in zip(tuple(interv_start_frame), tuple(interv_end_frame)):
        interv_index_bool.append(not all(f0_array[x:y] == 0))

# para quedarnos solo con las intervenciones que tienen valores de f0
    interv_start_frame_f0 = interv_start_frame[interv_index_bool]
    interv_end_frame_f0 = interv_end_frame[interv_index_bool]

    IOIstart_intervention_timeFrame = np.array(interv_start_frame_f0)
    IOIend_intervention_timeFrame = np.array(interv_end_frame_f0)

# inicio y final (s) de las intervenciones con f0 para una funcion posterior (get_interv_wav) 
# que segmenta el audio completo en intervenciones y las guarda 
    """START AND END TIME of voiced interventions"""
    interv_x_index = [time[i] for i in interv_start_frame_f0]
    interv_y_index = [time[i] for i in interv_end_frame_f0-1]

# Regresa los frames de las intervenciones    
    return IOIstart_intervention_timeFrame, IOIend_intervention_timeFrame, interv_x_index, interv_y_index


def intensity_peak_activity_frame(sound, min_pitch, samplerate, threshold):
    
    #############################
    #### ALL ABOUT INTENSITY ####   
    #############################
# utiliza la misma funcion que antes (to_intensity)
# y regresa el tiempo y la intensidad de todos los picos

    intensity = sound.to_intensity(int(min_pitch))
    intensity_peaks = call(intensity, "To IntensityTier (peaks)")
    intensity_peak_table = call(intensity_peaks, "Down to TableOfReal")
    
    npeaks = call(intensity_peak_table, "Get number of rows")
    peak_time = []
    peak_intensity = []
    for i in range(1, npeaks+1):
        peak_time.append(call(intensity_peak_table, "Get value", i, 1)) #calls the time of the intensity peaks
        peak_intensity.append(call(intensity_peak_table, "Get value", i, 2)) #calls the intensity

    peak_time_thres = []
    peak_intensity_thres = []
    for (x, y) in zip(tuple(peak_time), tuple(peak_intensity)):
        if y > threshold:
            peak_time_thres.append(x)
            peak_intensity_thres.append(y)

# se obtiene el frame de dichos picos
    IOI_intensity_timeFrame = np.array([int(i*samplerate) for i in peak_time_thres])
    
    return IOI_intensity_timeFrame

### Onset Detection
def onsetFunctions(sound, numOdfs, hopTime=0.005, int_threshold=-34):
    """Calculate onset detection functions with input parameters:
         sound (parselmouth.Sound): data to analyze
         numOdfs (int): The number of odfs to compute
         hopTime: hop size (in seconds) (time between successive frames)
       Returns value:
         odf: a matrix of onset detection function results, with one row for each."""
    
    odf = np.zeros((sound.values.shape[1], numOdfs))
    
    f_f0array = f0_array(sound,timeStep=hopTime)
    
    #IOI with intervention onsets
    f_interventionactivity = intervention_activity_frame(sound, f_f0array[0], f_f0array[1], f_f0array[3], int_threshold=int_threshold) #f0_array, min_pitch, time
    odf[f_interventionactivity[0],0] = 1
    
    #IOI with vowel onsets
    f_f0activity = f0_activity_frame(f_f0array[0], f_f0array[3]) #f0_array, time
    odf[f_f0activity[0],1] = 1
    
    return odf

def getOnsets(sound, numOfODFs=3, hop_time=0.005, int_threshold=-34):
    """Computes ODFs for a .wav file 
        sound (parselmouth.Sound): data to analyze
        numOfODFs (int): The number of ODFs to compute
        int_threshold (float): magnitude threshold in negative dB to consider a change in intervention
        
       Returns (array of onsets): the times at which onsets occur (at ms resolution. i.e. a value of 1 corresponds to 1 ms)"""
    
    odf = onsetFunctions(sound, numOfODFs, hopTime=hop_time, int_threshold=int_threshold)
    
    # standardise to zero mean, unit stdev
    odf = np.divide(np.subtract(odf, np.mean(odf)), np.std(odf))
    odf[np.isnan(odf)]=0
    
    t = np.multiply(range(len(odf)), 1/sound.get_sampling_frequency()) # <--- UNNECESARY?
    d = np.diff(odf, axis=0)
    isPeak = np.multiply(np.greater(np.concatenate([np.zeros((1,numOfODFs)), d]), 0),
                np.less(np.concatenate([d, np.zeros((1,numOfODFs))]), 0))
    peaks = np.nonzero(isPeak) #peaks is a tuple of arrays. First array is where the peaks are. Second is which odf.
    peakIndex = peaks[0]
    odfNum = peaks[1]
    
    return ((peakIndex/sound.get_sampling_frequency())*1000).astype(int) #return offsets in milliseconds

### Clustering algorithm using onsets/IOIs.
class Cluster:
    """Class to represent a cluster of IOIs."""
    def __init__(self, IOIs):
        """ Initiate with a list of IOIs"""
        self.IOIs = IOIs
        self.size = len(self.IOIs)
        self.interval = np.mean(self.IOIs)
        self.score = 0
    
    def add_IOIs(self, IOIs):
        """ Add IOIs to it's list and recomputes the length and interval."""
        self.IOIs = self.IOIs + IOIs
        self.size = len(self.IOIs)
        self.interval = np.mean(self.IOIs)

def scoreFor(d):
    """Cluster scoring function for integer (d) multiples of other clusters' intervals."""
    assert(type(d) == int)
    if d >= 1 and d <= 4:
        return 6-d
    elif d >= 5 and d <= 8:
        return 1
    else:
        return 0

def calculateClusters(onsets, clusterWidth=3, lowThreshold=100, highThreshold=1500):
    """ The clustering algorithm.
    1. Find all possible IOIs within a defined period (lowThreshold<IOI<highThreshold).
    2. Find clusters of all the IOIs.
    3. Give clusters scores based on their weight and how many other clusters divide into them.
    
    Returns: a list of tuples of the form: (Cluster interval, Cluster score)"""
    
    clusters = []

    for onset_1 in onsets:
        for onset_2 in onsets:
            IOI = onset_2 - onset_1
            added = False
            if IOI >= lowThreshold and IOI <= highThreshold:
                #IOIs should be greater than 100msec or less than 1500msec
                for cluster in clusters:
                    if added == False and np.abs(cluster.interval - IOI) < clusterWidth:
                        cluster.add_IOIs([IOI])
                        added = True
                if added == False:
                    clusters = clusters + [Cluster([IOI])]
    #Join clusters that have similar intervals.
    for cluster_1 in clusters:
        for cluster_2 in clusters:
            if cluster_1 != cluster_2:
                if np.abs(cluster_1.interval - cluster_2.interval) < clusterWidth:
                    cluster_1.add_IOIs(cluster_2.IOIs)
                    clusters.remove(cluster_2)
    
    #Score clusters based on its interval relative to other clusters' intervals.
    maxClusterScore = 0
    for cluster_1 in clusters:
        for cluster_2 in clusters:
            if cluster_1 != cluster_2:
                for n in range(8):
                    if np.abs(cluster_1.interval - n*cluster_2.interval) < clusterWidth:
                        cluster_1.score += scoreFor(n)*cluster_2.size
        if cluster_1.score > maxClusterScore:
            maxClusterScore = cluster_1.score
    
    #Normalise the cluster scores to between 0 and 1.
    if maxClusterScore > 0:
        for cluster in clusters:
            cluster.score = cluster.score/maxClusterScore
    
    #Return two lists: (cluster IOIs, cluster scores)
    clusterIntervals = []
    clusterScores = []
    for cluster in clusters:
        clusterIntervals = clusterIntervals + [cluster.interval]
        clusterScores = clusterScores + [cluster.score]
    
    #A sorted list of tuples of all the clusters - (IOI, score)
    sortedClusters = sorted(list(zip(clusterIntervals, clusterScores)), key=lambda tup: tup[0])
    
    return sortedClusters

### Multi-Agent Beat Tracking
class Agent():
    def __init__(self, _period, _onsetTime, _alpha, _clusterScore):
        self.period = _period
        self.onsetTime = _onsetTime
        self.alpha = _alpha
        self.mostRecentBeatTime = _onsetTime
        self.beatTimes = [_onsetTime]
        self.error = 0
        self.notFinished = True
        self.missedBeats = 0
        self.totalMissedBeats = 0
        self.hits = 0
        self.gotLost = False
        self.clusterScore = _clusterScore
        self.score = _clusterScore
    
    def predict(self):
        predictedBeat = self.mostRecentBeatTime + self.period
        return predictedBeat
    
    def updateHit(self, _nextBeat):
        #If the agent has become too fast (quicker than 150ms), then stop it.
        if _nextBeat - self.mostRecentBeatTime < 150:
            self.notFinished = False
            self.gotLost = True
        self.beatTimes = self.beatTimes + [_nextBeat]
        self.mostRecentBeatTime = _nextBeat
        self.missedBeats = 0
        self.hits = self.hits + 1
    
    def updateSemiHit(self, _nextBeat, _error):
        diff = _nextBeat - self.mostRecentBeatTime
        self.beatTimes = self.beatTimes + [_nextBeat]
        self.mostRecentBeatTime = _nextBeat
        self.period = self.alpha*diff + (1-self.alpha)*self.period
        self.error += (_error*_error)
        self.missedBeats = 0
        self.hits = self.hits + 1
    
    def updateMissedBeat(self):
        self.missedBeats = self.missedBeats+1
        self.totalMissedBeats = self.totalMissedBeats+1
        if self.missedBeats >= 4:
            self.notFinished = False
            self.gotLost = True
        predictedBeat = self.predict()
        self.mostRecentBeatTime = predictedBeat
        self.beatTimes = self.beatTimes + [predictedBeat]
    
    def getAveragePeriod(self):
        return ((self.beatTimes[-1:]-self.beatTimes[0])/len(self.beatTimes))[0]
        
    def saveBeats(self,textfilename):
        with open(textfilename, 'w') as filehandle:
            filehandle.writelines("%s\n" % str(int(beat)/1000) for beat in self.beatTimes)
    
    def calculateScore(self, scoreMatrix):
        if self.hits == 0:
            self.score = -999
            return
        if self.gotLost:
            self.score = -999
            return
        self.score = scoreMatrix[0]*np.log(1+self.hits)-scoreMatrix[1]*np.log(1+self.totalMissedBeats)+scoreMatrix[2]*self.clusterScore-scoreMatrix[3]*((self.error)/self.hits)
    

def findClosestOnset(uniqueOnsets, predictedBeatTime):
    #Subtract the predicted beat time. The minimum absolute value of the difference is the closest onset.
    difference = np.subtract(uniqueOnsets, predictedBeatTime)
    absDifference = np.abs(difference)
    return uniqueOnsets[np.argmin(absDifference)]


def letAgentsLoose(clusters, numberOfStartingOnsets, uniqueOnsets,
                   alpha, innerErrorThreshold, outerErrorThreshold):
    """Algorithm that creates a series of agents that track the beat over time.
        clusters (list of tuples: [(ClusterIOI, ClusterScore)]): The clusters that have been found.
        numberOfStartingOnsets (int): How many onsets we consider as possible start of an agent's beat
        uniqueOnsets: list of unique onsets found by onset detection
        alpha (0-1): determines an agent's momentum vs reactiveness to beats that don't match up exactly with it's prediction
        innerErrorThreshold: The window either side of an agents prediction that it will accept
                                an onset as part of the beat
        outerErrorThreshold: The outer window either side of an agents prediction where an onset will affect
                                the agent's period."""

    agents = []
    startingOnsets = uniqueOnsets[0:numberOfStartingOnsets]
    for i,onset in enumerate(startingOnsets):
        #Create an agent for each cluster for each onset.
        for cluster in clusters:
            period = cluster[0]
            clusterScore = cluster[1]
            agents = agents + [Agent(period, onset, alpha, clusterScore)]
    
    for t in range(150): #150 is enough time steps assuming a 300bpm maximum for 30s clips
        for agent in agents:
            if agent.notFinished:
                predictedBeatTime = agent.predict()
                closestOnset = findClosestOnset(uniqueOnsets, predictedBeatTime)
                #If the agent has reached the end then we don't want to carry on predicting for it
                if closestOnset == uniqueOnsets[len(uniqueOnsets)-1]:
                    agent.notFinished = False
                if agent.notFinished:
                    #Error between agent's prediction and closest onset
                    error = np.abs(closestOnset-predictedBeatTime)
                    if error < innerErrorThreshold:
                        #Move the agent forward one step.
                        agent.updateHit(closestOnset)
                    elif error < outerErrorThreshold:
                        #agent's period will update slightly based on this onset.
                        agent.updateSemiHit(closestOnset, error)
                    else:
                        #Charge forth anyways with a missed beat.
                        agent.updateMissedBeat()
    return agents

def getBestAgent(agents, scoreMatrix):
    """Look through all the agents and find the one that scored best.
        agents (list of agents)
        scoreMatrix ([a, b, c, d]): list of scoring weightings for each possible scoring items within an agents
                                    getScore() method."""
    bestAgent = -1
    bestScore = -99999999999999
    for agent in agents:
        agent.calculateScore(scoreMatrix)
        if (agent.score > bestScore):
            bestScore = agent.score
            bestAgent = agent
    return bestAgent

### Beatroot class
class Beatroot:
    """Class that runs the Beatroot algorithm."""
    def __init__(self,
                snd,
                rate,
                numOfODFs=2,
                clusterWidth=3,
                alpha=0.35,
                innerErrorThreshold=5,
                outerErrorThreshold=50,
                hop_time=0.005,
                int_threshold = -34,
                numberOfStartingOnsets=25,
                scoreMatrix = [10, 3, 5, 0.01]):
        """ Initiate Beatroot class:
        snd (array): data to analyze, read using soundfile.read
        rate (int): sample rate of data, read using soundfile.read
        clusterWidth (int): size of clusters
        alpha (float, 0-1): determines an agent's momentum vs reactiveness to beats that don't match up exactly with it's prediction
        innerErrorThreshold (int): The window in ms either side of an agents prediction that it will accept an onset as part of the beat
        outerErrorThreshold (int): The outer window in ms either side of an agents prediction where an onset will affect the agent's period.
        hop_time (float): minimum hop in time (s) between peaks.
        int_threshold (float): magnitude threshold in negative dB to consider a change in intervention
        numberOfStartingOnsets (int): Number of agents is numberOfStartingOnsets x numberOfClusters.
        scoreMatrix ([gHitScore, gMissScore, gSalienceScore, gErrorScore]): variables that determine an agent's score.
        """
        self.snd = snd
        self.rate = rate
        self.numOfODFs = numOfODFs
        self.clusterWidth = clusterWidth
        self.alpha = alpha
        self.innerErrorThreshold = innerErrorThreshold
        self.outerErrorThreshold = outerErrorThreshold
        self.hop_time = hop_time
        self.int_threshold=int_threshold
        self.numberOfStartingOnsets = numberOfStartingOnsets
        self.scoreMatrix = scoreMatrix
    
    def run(self):
        """Test beat detection, printing useful diagnostics."""
        sound = parselmouth.Sound(self.snd, self.rate)
        sound = sound.resample(new_frequency = 16000)
        
        #get onsets from praat algorithms
        onsets = getOnsets(sound, self.numOfODFs, self.hop_time, self.int_threshold)
        
        #initialize clusters
        clusters = calculateClusters(onsets, self.clusterWidth)
        
        #Find the ordered set of onset times. Necessary when multiple ODFs used.
        uniqueOnsets = sorted(set(onsets))
        
        #run multi-agent search algorithm to find rhythm
        agents = letAgentsLoose(clusters, self.numberOfStartingOnsets, uniqueOnsets, self.alpha, self.innerErrorThreshold, self.outerErrorThreshold)
        
        #calculate the best agent
        bestAgent = getBestAgent(agents.copy(), self.scoreMatrix)
        
        if bestAgent == -1:
            return -1
        else:
            #return bestAgent.getAveragePeriod()
            return 60/(bestAgent.period/1000) #returning in the form of BPM


if __name__ == "__main__":
    #argument: wav file
    
    wave_files = 'EN2002b.Mix-Headset.wav'
    
    if len(sys.argv) > 1:
        wave_files = sys.argv[1]
    #else:
    #    exit()
    
    rate = 16000
    minimum_intervention_length_secs = 2.5
    int_threshold = -34
    timeStep=0.5
    minpause=2.0
    
    json_basename=("".join(x for x in os.path.basename(wave_files) if x.isalnum() or x in "._-")).replace('.wav','')
    
    audio_onset = {}
    wave_list = glob.glob(wave_files)
    for wave_file in wave_list:
        print("Analyzing",wave_file,":")
        if not os.path.isfile(wave_file):
            print("File '", wave_file,"' does not exists. Exiting.")
            break
        
        print("Gathering info from WAV file to initialize variables...")
        file_snd, file_rate = soundfile.read(os.path.join(wave_file))
        
        print("Extracting BPM from interventions...")
        wave_file_basename = os.path.basename(wave_file)
        #print(wave_file_basename)
        
        #obtaining the set of interventions
        sound = parselmouth.Sound(file_snd,file_rate)
        sound = sound.resample(new_frequency = rate)
        f_f0array = f0_array(sound,timeStep=timeStep)
        
        #IOI with intervention onsets
        f_interventionactivity = intervention_activity_frame(sound, f_f0array[0], f_f0array[1], f_f0array[3], int_threshold = int_threshold, minpause=minpause) #f0_array, min_pitch, time
        #print(f_interventionactivity)
        
        #storing onset information for this wav file
        audio_onset[wave_file_basename] = {}
        audio_onset[wave_file_basename]["start"] = f_interventionactivity[0].tolist()
        audio_onset[wave_file_basename]["end"] = f_interventionactivity[1].tolist()
        audio_onset[wave_file_basename]["bpm"] = [-1] * len(audio_onset[wave_file_basename]["start"])
        
        for onset_i, onset_start in enumerate(audio_onset[wave_file_basename]["start"]):
            onset_end = audio_onset[wave_file_basename]["end"][onset_i]
            
            if (onset_end-onset_start)/rate > minimum_intervention_length_secs:
                snd = file_snd[onset_start:onset_end]
                beatroot = Beatroot(snd,rate)
                try:
                    this_bpm = beatroot.run()
                    if this_bpm > 0 and this_bpm < 200:
                        audio_onset[wave_file_basename]["bpm"][onset_i] = this_bpm
                        print("    "+str(onset_start)+" to "+str(onset_end)+" -> "+str(audio_onset[wave_file_basename]["bpm"][onset_i]))
                    else:
                        print("    "+str(onset_start)+" to "+str(onset_end)+" -> no good agent")
                        audio_onset[wave_file_basename]["bpm"][onset_i] = -1
                        audio_onset[wave_file_basename]["start"][onset_i] = -1
                        audio_onset[wave_file_basename]["end"][onset_i] = -1
                except:
                    print("    "+str(onset_start)+" to "+str(onset_end)+" -> could not beatroot")
                    audio_onset[wave_file_basename]["bpm"][onset_i] = -1
                    audio_onset[wave_file_basename]["start"][onset_i] = -1
                    audio_onset[wave_file_basename]["end"][onset_i] = -1
            else:
                print("    "+str(onset_start)+" to "+str(onset_end)+" -> too short")
                audio_onset[wave_file_basename]["bpm"][onset_i] = -1
                audio_onset[wave_file_basename]["start"][onset_i] = -1
                audio_onset[wave_file_basename]["end"][onset_i] = -1
        
        
        audio_onset[wave_file_basename]["bpm"] = [value for value in audio_onset[wave_file_basename]["bpm"] if value != -1]
        audio_onset[wave_file_basename]["start"] = [value for value in audio_onset[wave_file_basename]["start"] if value != -1]
        audio_onset[wave_file_basename]["end"] = [value for value in audio_onset[wave_file_basename]["end"] if value != -1]
        
        onset_file = open(json_basename+"__bpm.json", "w")
        json.dump(audio_onset, onset_file)
        onset_file.close()
        
