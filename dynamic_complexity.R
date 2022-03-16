####################### Introduction ###########################################
# Article: Structural complexity of the long-term collective rhythm of 
#          naturalistic conversations.
# Year: 2022
# Author: Arodi Farrera, Caleb Rascon, Gabriel Ramos-Fernández



########################## Summary #############################################
# Source code to generate the analyzes shown in Figure 1.A, 1.B, 1.C, 1.D, using
# one meeting as an example.
# This code is divided into three main sections
# 1. DATA: 
# 1.1. Acoustic Data: reading TXT files containing acoustic onsets for one 
#      meeting. These files come from THE PYTHON CODE
# 1.1.1. Types of time series: steps to generate the three types of time series:
#        Observed, Random, and Isochronous (see Figure 1.A).
# 1.2. BPM data: reading the TXT file containing the BPM values for intervention

# 2. METRICS:
# 2.1. Burstiness: code to compute the B parameter
# 2.2. RQA: code to compute %DET, %LAM, ENTR
# 2.3. Allan factor analysis: code to compute the exponent alpha

# 3. PLOTS:
# 3.1. 3D Plot: code to generate the 3D plot of median values for all 
#               types of time series and all meetings (Figure 2)
# 3.2. Over time: code to generate Figure 4 showing the dynamic structural
#                 organization and tempo characteristic of the selected meeting
# 3.3. Boxplots: code to generate Figure 3, but instead of comparing different
#                meetings, it compares the selected meeting against the
#                random and isochronous versions when appropriate.

## The functions used in section 2.METRICS were adapted from:
# Xu, L., de Barbaro, Abney, D., Cox, R. (2020). Finding Structure in Time: 
#     Visualizing and Analyzing Behavioral Time Series. Technology and Code. 
#     doi: 10.3389/fpsyg.2020.01457





## packages needed
library("dplyr")
library("tidyr")
library("pracma")
library("ggdist")
library("ggplot2")
library("crqa")
library("slider")
library("plotly")
library("egg")
library("patchwork")


#### 1.DATA ####
#### 1.1.Acoustic data ####
filenames <- dir(pattern = "*Mix-Headset.txt")
samplerate <- 16000
IOIstart <- read.table("FRAME_longintS_EN2002b.Mix-Headset.txt")
IOIend <- read.table("FRAME_longintE_EN2002b.Mix-Headset.txt")
IOIstartS <- read.table("FRAME_shortintS_EN2002b.Mix-Headset.txt")
IOIvowel <- read.table("FRAME_shortvowel_28589568_EN2002b.Mix-Headset.txt")

#### 1.1.1.Types of time series ####
set.seed(1)

## OBSERVED 
ioiDur <- mean(diff(IOIvowel$V1)/samplerate) #mean IOI duration
startOfPattern <- min(IOIvowel$V1)/samplerate #starting point for the random and periodic version
patIOI <- IOIvowel$V1/samplerate #moment in time (s.) with acoustic activity
nlength <- length(IOIvowel$V1)

## ISOCHRONOUS VERSION
xIso <- array(ioiDur, dim = nlength) #array of same length as the original for periodic version
# creating an isochronous pattern with IOIs = ioiDur, starting at startOfPattern
ioiSums <- cumsum(xIso) #cumulative sum
patIso <- c(startOfPattern, ioiSums[-nlength]+ startOfPattern) #isochronous version
# isochronous version with noise
patIsonoise <- c(patIso[1], patIso[-1] + rnorm(nlength-1, sd = 0.1))


### Time series with a window size of 1 s.
## OBSERVED
patron <- patIOI
time <- data.frame(time = seq(1, length(patron)))
# counts the number of onsets per second
seconds <- data.frame(table(floor(patron)))
colnames(seconds) <- c("time", "Sec")
seconds$time <- as.numeric(as.character(seconds$time))
patronLong <- max(seconds$time)
# Data frame in which unvoiced events = 0, and voiced events > 0
discIOI <- left_join(time,seconds, by = "time") 
discIOI <- discIOI[1:patronLong,]
discIOI[is.na(discIOI)] = 0 

## RANDOM VERSION
discIOIshuff <- discIOI
# the observed time series but randomly reordered 
discIOIshuff$Sec <- sample(discIOI$Sec)

## ISOCHRONOUS VERSION
#same process as for the observed time series but for the Isochronous version
patron <- patIsonoise
time <- data.frame(time = seq(1, length(patron)))
seconds <- data.frame(table(floor(patron)))
colnames(seconds) <- c("time", "Sec")
seconds$time <- as.numeric(as.character(seconds$time))
discIsonoise <- left_join(time,seconds, by = "time")
discIsonoise <- discIsonoise[1:patronLong,]
discIsonoise[is.na(discIsonoise)] = 0


#### 1.2.BPM data ####
BPM2b <- read.table("BPM_EN2002b.Mix-Headset.txt", header = T)
BPM2b$start.s <- BPM2b$start/samplerate # start in seconds
BPM2b$end.s <- BPM2b$end/samplerate # end in seconds

## Adding acceleration
BPM2b$acc <- c(NA, with(BPM2b, diff(bpm)/diff(start/samplerate)))
BPM2b$meeting <- "EN2002b"



#### 2.METRICS ####
#### 2.1.Burstiness ####
bursti <- function(subArray){
  # see the matlab code in XU et al. 2020
  bursti <- (sd(subArray)-mean(subArray))/(sd(subArray)+mean(subArray))
  return(bursti)
}

# Arguments for the sliding window
window_size = 60*2 #two minutes
step_size = 1 # one second

# the sliding window technique cannot be applied to the last 2 minutes
# because the step size is of one second, thus, we need 2 minutes of NA values
# at the end of the array of B estimates
NAs <- array(NA, dim = 60*2)

####  B metric
## OBSERVED
discIOI$B <- c(unlist(slide(discIOI$Sec, bursti, .after = window_size, 
      .step = step_size, .complete = T)), NAs)

## RANDOM VERSION
discIOIshuff$B <- c(unlist(slide(discIOIshuff$Sec, bursti, .after = window_size, 
                                 .step = step_size, .complete = T)), NAs)
## ISOCHRONOUS VERSION
discIsonoise$B <- c(unlist(slide(discIsonoise$Sec, bursti, .after = window_size, 
                                 .step = step_size, .complete = T)), NAs)


#### 2.2.Allan factor analysis ####
allanplotter <- function(fin, powers){
  # Adapted from the matlab scripts in Xu et al. (2020)
  # range of seconds used for each window size
  powers = powers
  base = 2   
  start = 2
  
  potencia = seq(start, powers, by = 0.5)
  numint = base^potencia  # Number of adjacent windows
  
  #ALLAN FACTOR  
  AllanF <- list()
  for (j in 1:length(numint)) {
    # I. get mean event count
    medio = sum(fin)/numint[j] # average of uniforminly distributed events given the number of adjacent windows
    
    #II. get Allan variance of spike count
    int = floor(length(fin)/numint[j]) # get length of each interval
    indez <- c(seq(1, int*numint[j], by = int), int*numint[j]) # index of where does each interval start
    
    # get variance
    tot <- c()
    for (k in 1:(length(indez)-2)) {
      # the differences in counts of events between adjacent windows squared
      tot[k] <- (sum(fin[ (indez[k]:indez[(k+1)])] ) - sum(fin[ indez[(k+1)] : indez[k+2] ] ))^2
    }
    AllanF[[j]] <- (sum(tot)/length(tot))/(2*medio) 
  }
  
  AllanF_outcome <- data.frame(AF = rev(unlist(AllanF)), 
                               abcisa = rev((length(fin)/numint))) #in seconds
  return(AllanF_outcome)
  
}
fractal <- function(data){
  powers = 8
  fractal_outcome <- allanplotter(as.array(data), powers)
  abcisa <- fractal_outcome$abcisa
  AT <- fractal_outcome$AF
  p_all <- polyfit(log(abcisa), log(AT), 1)[1]
  
  return(p_all)
} # for sliding window

# Arguments for the sliding window
window_size = 60*10 # a ten minutes window size to avoid a biased estimate
NAs <- rep(NA, times = 60*10) # same reason as before, but now with a 10 min. vector of NAs

####  Exponent alpha 
## OBSERVED
alpha <- unlist(as.array(slide(discIOI$Sec, fractal, .after = window_size, 
                               .step = step_size, .complete = T)))
discIOI$alpha <- c(alpha, NAs)

## RANDOM VERSION
alpha <- unlist(as.array(slide(discIOIshuff$Sec, fractal, .after = window_size, 
                               .step = step_size, .complete = T)))
discIOIshuff$alpha <- c(alpha, NAs)

## ISOCHRONOUS VERSION
alpha <- unlist(as.array(slide(discIsonoise$Sec, fractal, .after = window_size, 
                               .step = step_size, .complete = T)))
discIsonoise$alpha <- c(alpha, NAs)



#### 2.3.RQA ####
crqa_categorical <- function(data){
  #see Xu et al. (2020) for an introduction to the topic
  crqa_output <- crqa(as.array(data), as.array(data), datatype = "categorical")[c(1,2,6,8)]
  
  output <- list(ENTR = crqa_output$ENTR, # Entropy 
                 DET = crqa_output$DET,   # %DET
                 LAM = crqa_output$LAM)   # %LAM
  return(output) 
}


#### RQA measurentes 
## OBSERVED
recur <- bind_rows(slide(discIOI$Sec, crqa_categorical, .after = window_size, 
                         .step = step_size, .complete = T))
discIOI$entropy <- c(recur$ENTR, NAs) #adding NAs at the end as before
discIOI$DET  <- c(recur$DET, NAs)
discIOI$LAM <- c(recur$LAM, NAs)

## RANDOM VERSION
recur <- bind_rows(slide(discIOIshuff$Sec, crqa_categorical, .after = window_size, 
                         .step = step_size, .complete = T))
discIOIshuff$entropy <- c(recur$ENTR, NAs)
discIOIshuff$DET  <- c(recur$DET, NAs)
discIOIshuff$LAM <- c(recur$LAM, NAs)

## ISOCHRONOUS VERSION
recur <- bind_rows(slide(discIsonoise$Sec, crqa_categorical, .after = window_size, 
                         .step = step_size, .complete = T))
discIsonoise$entropy <- c(recur$ENTR, NAs)
discIsonoise$DET  <- c(recur$DET, NAs)
discIsonoise$LAM <- c(recur$LAM, NAs)







## Joining discIOI and BPM2b
# patron <- BPM2b$start.s
# time <- data.frame(time = seq(1, length(discIOI$time)))
# seconds <- data.frame(table(floor(patron)))
# seconds$BPM <- BPM2b$bpm 
# colnames(seconds) <- c("time", "Sec", "BPM")
# seconds$time <- as.numeric(as.character(seconds$time))
# timeBPM <- left_join(time,seconds, by = "time") 
# discIOI <- left_join(discIOI, timeBPM[,-2], by = "time")

obsEN2002b <- discIOI
obsEN2002b$meeting <- "EN2002b" 
obsEN2002b$type <- "Observed" 

ranEN2002b <- discIOIshuff 
ranEN2002b$meeting <- "EN2002b" 
ranEN2002b$type <- "Random" 

isoEN2002b <- discIsonoise
isoEN2002b$meeting <- "EN2002b" 
isoEN2002b$type <- "Isochronous" 

## Median values for observed, random and isochronous patterns, respectively
rbind(apply(obsEN2002b[,c(3,4,7)], 2, function(x) median(x, na.rm = T)),
 apply(ranEN2002b[,c(3,4,7)], 2, function(x) median(x, na.rm = T)),
 apply(isoEN2002b[,c(3,4,7)], 2, function(x) median(x, na.rm = T)))

## SD values for observed, random and isochronous patterns, respectively
rbind(apply(obsEN2002b[,c(3,4,7)], 2, function(x) sd(x, na.rm = T)),
apply(ranEN2002b[,c(3,4,7)], 2, function(x) sd(x, na.rm = T)),
apply(isoEN2002b[,c(3,4,7)], 2, function(x) sd(x, na.rm = T)))

#### 3.PLOTS ####
#### 3.1.3D PLOT ####
# Median values for all meetings
DF3D <- read.table("DF3D.txt", header = T)
DF3D$Type <- factor(DF3D$Type, levels = c("Observed", "Random", "Isochronous"))
## format used for displaying meeting names on the 3D plot
DF3D$names <- c("EN2002b", "", "", 
                "EN2001d", "EN2001e", "EN2001a", 
                "","", "", "", "", "", "","", "", "", "", "")

# colors
paleta <- hcl.colors(5, "Viridis")[c(2,4,5)]
f <- plot_ly(DF3D, 
             x = ~B, y = ~Entropy, z = ~Alpha, color = ~Type,
             marker = list(size = 10), 
             colors = paleta, text = DF3D$names)

f <- f %>% add_markers() %>% layout(
  scene = list(aspectratio = list(x = 1,y = 1,z = 1),
               xaxis = list(range = c(-1,1)),
               zaxis = list(range = c(-1,1))))

t = list(size = 12, color = toRGB("black"))
f <- f %>% add_text(textfont = t, textposition = "top right")

f %>% layout(scene = list(annotations = list(
  list(x = -0.16,y = 0.5,z = 0.37,
       showarrow = F,
       text = "EN2001b",
       xanchor = "top left",
       font = list(color = "black",size = 12)), 
  list(x = -0.1,y = 0.62,z = 0.5,
       showarrow = F,
       text = "EN2002a",
       xanchor = "top left",
       font = list(color = "black",size = 12))) ))



#### 3.2.Boxplots ####
# B metric compared to random and periodic versions
e <- ggplot(DFtime, 
            aes(x = type, y = B, fill = type))+
  geom_boxplot() + theme_bw() + ylim(c(-1, 1)) +
  theme(strip.background = element_blank(),
        strip.text.y = element_blank(),
        axis.title.x=element_blank(),
        legend.position = "none") + ylab("B") +
  scale_fill_manual(values=c(paleta))

# Exponent alpha compared to random and periodic versions
f <- ggplot(DFtime, aes(x = type, y = alpha, fill = type))+
  geom_boxplot() + theme_bw() + ylim(c(-1, 1)) +
  theme(strip.background = element_blank(),
        strip.text.y = element_blank(),
        axis.title.x=element_blank(),
        legend.position = "none") + ylab("Alpha") +
  scale_fill_manual(values=c(paleta))

# Entropy compared to random and periodic versions
g <- ggplot(DFtime, aes(x = type, y = entropy, fill = type))+
  geom_boxplot() + theme_bw() + ylim(c(0, 5)) +
  theme(strip.background = element_blank(),
        strip.text.y = element_blank(),
        axis.title.x=element_blank(),
        legend.position = "none") + ylab("Entropy") +
  scale_fill_manual(values=c(paleta))

# %DET compared to random and periodic versions
h <- ggplot(DFtime, aes(x = type, y = DET, fill = type))+
  geom_boxplot() + theme_bw() + ylim(c(0, 100)) +
  theme(strip.background = element_blank(),
        strip.text.y = element_blank(),
        axis.title.x=element_blank(),
        legend.position = "none") + ylab("%DET") +
  scale_fill_manual(values=c(paleta))

# %LAM compared to random and periodic versions
i <- ggplot(DFtime, aes(x = type, y = LAM, fill = type))+
  geom_boxplot() + theme_bw() + ylim(c(0, 100)) +
  theme(strip.background = element_blank(),
        strip.text.y = element_blank(),
        axis.title.x=element_blank(),
        legend.position = "none") + ylab("%LAM") +
  scale_fill_manual(values=c(paleta))

# BPM 
j <- ggplot(BPM2b, aes(x = meeting, y = bpm, fill = meeting))+
  geom_boxplot() + theme_bw() + ylim(c(40, 200)) +
  theme(strip.background = element_blank(),
        strip.text.y = element_blank(),
        axis.title.x=element_blank(),
        legend.position = "none") + ylab("BPM") +
  scale_fill_manual(values=c(paleta[1]))  

etoj <- (e + f + g ) / (h + i + j) + plot_layout(ncol = 1, guides = "collect")
# Saving the plot
file_name = paste("./Summary", ".png", sep = "")
etoj
ggsave(file_name, device = "tiff", width=15, height=10, dpi=150, units="in", compression = "lzw")


#### 3.3.Over time: b, alpha, entropy, acceleration ####
head(obsEN2002b)
DFtime <- rbind(obsEN2002b, ranEN2002b, isoEN2002b)
DFtime$type <- factor(DFtime$type, levels = c("Observed", "Random", "Isochronous"))

paleta <- hcl.colors(5, "Viridis")[c(2, 4, 5)]
# B metric over time
a <- ggplot(DFtime, aes(x = time, y = B, color = type)) +
  geom_line(size = 0.5) + ylim(c(-1, 1)) + xlim(c(0, 1700))+
  theme_bw() +
  theme(strip.background = element_blank(),
        strip.text.y = element_blank(),
        axis.title.x=element_blank(),
        text = element_text(size = 30)) + 
  scale_colour_manual(values = paleta) + ylab("B")

a <- tag_facet(a, 
               x = -Inf, y = 1, 
               open = "", close = "",
               size = 10,
               tag_pool = "EN2002b")

# Exponent alpha over time
b <- ggplot(DFtime, aes(x = time, y = alpha, color = type)) +
  geom_line(size = 0.5) + ylim(c(-1, 1)) + xlim(c(0, 1700))+
  theme_bw() +
  theme(strip.background = element_blank(),
        strip.text.y = element_blank(),
        axis.title.x=element_blank(),
        text = element_text(size = 30)) + 
  scale_colour_manual(values = paleta) + ylab("Alpha")

# Entropy over time
c <- ggplot(DFtime, aes(x = time, y = entropy, color = type)) +
  geom_line(size = 0.5) + ylim(c(0, 6)) + xlim(c(0, 1700))+
  theme_bw() +
  theme(strip.background = element_blank(),
        strip.text.y = element_blank(),
        axis.title.x=element_blank(),
        text = element_text(size = 30)) + 
  scale_colour_manual(values = paleta) + ylab("Entropy")

# BPM acceleration over time
d <- ggplot(BPM2b, aes(x = start/samplerate, y = acc)) +
  geom_line(size = 0.5, color = paleta[1]) + xlim(c(0, 1700))+
  ylim(c(-20, 25))+ 
  theme_bw() +
  theme(strip.background = element_blank(),
        strip.text.y = element_blank(),
        text = element_text(size = 30)) + 
  scale_colour_manual(values = paleta[1]) +
  xlab("Time (s.)") + ylab("Acceleration")



abcd <- (a + b + c + d) + plot_layout(ncol = 1, guides = "collect")
# Saving the plot
file_name = paste("./Overtime", ".png", sep = "")
abcd
ggsave(file_name, device = "tiff", width=15, height=20, dpi=150, units="in", compression = "lzw")


