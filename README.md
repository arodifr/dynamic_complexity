# DYNAMIC COMPLEXITY
Source code to generate all the analyzes of the article "Structural complexity of the long-term collective rhythm of naturalistic conversations"


Applying the following workflow to the audios of the six analyzed meetings (EN2001d, EN2001e, EN2001a, EN2002b, EN2002a and EN2001b) will deliver the results presented in the article. These audios are part of the AMI Corpus (Carletta et al., 2006), which can be found at https://groups.inf.ed.ac.uk/ami/corpus/.   
   

![This is an image](https://github.com/arodifr/dynamic_complexity/blob/main/Figure1.png)

To make the code more readable, we exemplify this workflow using only one of these meetings (EN2002b). This workflow requires using Python version 3.8.8 and R studio version 4.3.0. 

### 1. **acoustic_onsets.py**    
Source code to generate the TXT files containing acoustic onsets (Figure 1.A) used as inputs in dynamic_complexity.R    
- Input: EN2002b.Mix-Headset.WAV (audio file to be analyzed)
- Output: FRAME_shortvowel.txt (f0 activity onsets)


### 2. **Beatroot.py**    
Source code to generate the TXT files containing BPM values (Figure 1.E, 1.F) used as inputs in dynamic_complexity.R
- Input: EN2002b.Mix-Headset.WAV (audio file to be analyzed)
- Output:	BPM.txt (start and end of interventions where BPM was estimated, and BPM estimates)


### 3. **dynamic_complexity.R**    
Source code to generate the analyzes shown in Figure 1.A, 1.B, 1.C, 1.D, 1.G, and 1.H using one meeting as an example.    
- Inputs:    
•	FRAME_shortvowel.txt (acoustic_onsets.py output)    
•	BPM.txt (Beatroot.py output)    
- Outputs:    
•	Summary.png: Boxplot of all metrics (B, exponent alpha, RQA measurements, BPM) obtained for the meeting used as an example and its random and isochronous versions.    
•	Overtime.png: Time series of four metrics (B, exponent alpha, entropy, and BPM’s acceleration) obtained from the meetings used as example.   

### 4. **dynamic_complexity[2nd].R**    
Source code to generate the analyzes shown in Figure 1.A, 1.B, 1.C, 1.D, 1.G, and 1.H using all meetings.    
- Inputs:    
•	All_metrics.txt (dynamic_complexity.R outputs, using all meetings)    
•	All_BPM.txt (Beatroot.py outputs, using all meetings)    
- Outputs:    
•	Figure2.png: 3D plot of median values of three metrics: B, exponent alpha, entropy.    
•	Figure3.png: Boxplot of all the metrics analyzed.    
•	Figure4.png: Loess smooth curve with 95% confidence intervals depicting the evolution of the structural complexity.    

### Additional files:    
- Figure1: schematic representation of the methodology.    

#### Reference
Carletta J, Ashby S, Bourban S, Flynn M, Guillemot M, Hain T, et al. The AMI Meeting Corpus: A Pre-announcement. In: Renals S, Bengio S, editors. Machine Learning for Multimodal Interaction. Berlin, Heidelberg: Springer Berlin Heidelberg; 2006.

