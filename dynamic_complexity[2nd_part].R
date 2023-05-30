####################### Introduction ###########################################
# Article: Structural complexity of the long-term collective rhythm of 
#          naturalistic conversations.
# Year: 2023
# Author: Arodi Farrera, Caleb Rascon, Gabriel Ramos-Fernandez



########################## Summary #############################################
# Source code to generate Figure 2, 3, and 4
# This code is divided into three main sections
# 1. DATA: 
# 1.1.Burstiness, Entropy, DET, LAM, Alpha: reading TXT file containing 
#      these metrics for each meeting 
#      This file comes from dynamic_complexity_2023.R
# 1.2. BPM data: reading the TXT file containing the BPM and acceleration values 
#      for each dialogue intervention within the meetings.

# 2. PLOTS:
# 2.1. 3D plot (Figure 2): code to generate the 3D plot using the median values 
#      for all types of time series and all meetings. 
# 2.2. Boxplots (Figure 3): code to generate a plot showing the distribution
#      of each of the metrics for each metric
# 2.3. LOESS (Figure 4): code to compare the trend over time of the structural
#      organization in Type I and Type II groups



## packages needed
library(plotly)
library(dplyr)
library(ggplot2)
library(patchwork)

#### 1.DATA ####
#### 1.1.Burstiness, Entropy, DET, LAM, Alpha ####
structureDF <- read.table("All_metrics.txt", header = T)
colnames(structureDF) <- c("Time", "B", "Entropy", "DET", "LAM", "Alpha", "Meeting", "Type")

#### 1.2.BPM ####
bpmDF <- read.table("All_BPM.txt", header = T)
bpmDF$Meeting <- factor(bpmDF$Meeting, levels = c("EN2002b", "EN2002a", "EN2001b",
                                                  "EN2001d", "EN2001e", "EN2001a"))
bpmDF$Time <- floor(bpmDF$start) 

# Labeling Type I and Type II groups
bpmDF <-bpmDF %>%
  mutate(Type = ifelse(Meeting == "EN2002b" | Meeting == "EN2002a" | Meeting == "EN2001b", "Type I", "Type II"))

#### 1.3. Observed time series ####
observedDF <- structureDF[structureDF$Type == "Observed" , ]
observedDF$Meeting <- factor(observedDF$Meeting, levels=c("EN2002b", "EN2002a", "EN2001b",
                                                              "EN2001d", "EN2001e", "EN2001a"))
# Labeling Type I and Type II groups
observedDF <-observedDF %>%
  mutate(Type = ifelse(Meeting == "EN2002b" | Meeting == "EN2002a" | Meeting == "EN2001b", "Type I", "Type II"))


#### 2.PLOTS ####
##### 2.1. 3D ####
# Median values
data3D <- structureDF[,c("B","Entropy","Alpha", "Meeting", "Type")] %>%
  group_by(Meeting, Type) %>%
  summarise(meanB = median(B, na.rm = T),
            meanEntropy = median(Entropy, na.rm = T),
            meanAlpha = median(Alpha, na.rm = T))
colnames(data3D) <- c("Meeting", "Type", "B", "Entropy", "Alpha")

######  Output: Table 1supp ####
data3D

# Plot
data3D$Type <- factor(data3D$Type, levels = c("Observed", "Random", "Isochronous"))
## format used for displaying meeting names on the 3D plot
data3D$names <- c("", "EN2001a", "", 
                  "", "EN2001b", "", 
                  "","EN2001d", "", 
                  "", "EN2001e", "", 
                  "","EN2002a", "", 
                  "", "EN2002b", "")
# colors
paleta <- hcl.colors(5, "Viridis")[c(2,4,5)]
f <- plot_ly(data3D, 
             x = ~B, y = ~Entropy, z = ~Alpha, color = ~Type,
             marker = list(size = 10), 
             colors = paleta, text = data3D$names
             )

f <- f %>% add_markers() %>% layout(
  scene = list(aspectratio = list(x = 1,y = 1,z = 1),
               xaxis = list(range = c(-1,1)),
               zaxis = list(range = c(-1,1))))

t = list(size = 12, color = toRGB("black"))
f <- f %>% add_text(textfont = t, textposition = "top right")

######  Output: Figure 2 ####
f

##### 2.2. Boxplots ####
# Burstiness
a <- ggplot(observedDF, aes(x = Meeting, y = B))+
  geom_boxplot(fill=paleta[1]) + 
  theme_bw() + ylim(c(-1, 1)) +
  theme(strip.background = element_blank(),
        strip.text.y = element_blank(),
        axis.title.x=element_blank(),
        legend.position = "none") + ylab("B") 

# Exponent alpha
b <- ggplot(observedDF, aes(x = Meeting, y = Alpha))+
  geom_boxplot(fill=paleta[1]) + theme_bw() + ylim(c(-1, 1)) +
  theme(strip.background = element_blank(),
        strip.text.y = element_blank(),
        axis.title.x=element_blank(),
        legend.position = "none") + ylab("Alpha") 

# Entropy
c <- ggplot(observedDF, aes(x = Meeting, y = Entropy))+
  geom_boxplot(fill=paleta[1]) + theme_bw() + ylim(c(0, 5)) +
  theme(strip.background = element_blank(),
        strip.text.y = element_blank(),
        axis.title.x=element_blank(),
        legend.position = "none") + ylab("Entropy") 

# %DET
d <- ggplot(observedDF, aes(x = Meeting, y = DET))+
  geom_boxplot(fill=paleta[1]) + theme_bw() + ylim(c(0, 100)) +
  theme(strip.background = element_blank(),
        strip.text.y = element_blank(),
        axis.title.x=element_blank(),
        legend.position = "none") + ylab("%DET") 

# %LAM 
e <- ggplot(observedDF, aes(x = Meeting, y = LAM))+
  geom_boxplot(fill=paleta[1]) + theme_bw() + ylim(c(0, 100)) +
  theme(strip.background = element_blank(),
        strip.text.y = element_blank(),
        axis.title.x=element_blank(),
        legend.position = "none") + ylab("%LAM") 

# BPM 
f <- ggplot(bpmDF, aes(x = Meeting, y = BPM))+
  geom_boxplot(fill=paleta[1]) + theme_bw() + ylim(c(40, 200)) +
  theme(strip.background = element_blank(),
        strip.text.y = element_blank(),
        axis.title.x=element_blank(),
        legend.position = "none") + ylab("BPM") 

atof <- (a + b + c ) / (d + e + f) + plot_layout(ncol = 1, guides = "collect")
######  Output: Figure 3 ####
atof

file_name = paste("./Meetings", ".png", sep = "")
atof
ggsave(file_name, device = "tiff",
       width=20, height=15, dpi=150, units="in", compression = "lzw")

##### 2.3. LOESS ####
# Burstiness
a <- na.omit(observedDF)  %>%
  ggplot(aes(x = Time, y = B, color = Type, group = Meeting)) +
#  geom_point(alpha = 0.05)+
  geom_smooth(aes(color = Type, fill = Type), method = "loess", se = T,
              alpha = .45) +
  theme_bw() +
  ylim(c(-1, 1)) +
  xlim(0, 5100) +
  theme(strip.background = element_blank(),
        strip.text.y = element_blank(),
        axis.title.x = element_blank(),
        text = element_text(size = 20)) + 
  theme(legend.key.size = unit(1.5, 'cm'),
        legend.text = element_text(size=20))+
  ylab("B")+
  scale_color_manual(values = c("#808000", "#FFC200")) +
  scale_fill_manual(values = c("#808000", "#FFC200"))

# Alpha
b <- na.omit(observedDF)  %>%
  ggplot(aes(x = Time, y = Alpha, color = Type, group = Meeting)) +
#  geom_point(alpha = 0.05)+
  geom_smooth(aes(color = Type, fill = Type), method = "loess", se = T,
              alpha = .45) +
  theme_bw() +
  ylim(c(-1, 1)) +
  xlim(0, 5100) +
  theme(strip.background = element_blank(),
        strip.text.y = element_blank(),
        axis.title.x = element_blank(),
        text = element_text(size = 20)) + 
  theme(legend.key.size = unit(1.5, 'cm'),
        legend.text = element_text(size=20))+ 
  ylab("Alpha")+
  scale_color_manual(values = c("#808000", "#FFC200")) +
  scale_fill_manual(values = c("#808000", "#FFC200"))

# Entropy
c <- na.omit(observedDF)  %>%
  ggplot(aes(x = Time, y = Entropy, color = Type, group = Meeting)) +
#  geom_point(alpha = 0.05)+
  geom_smooth(aes(color = Type, fill = Type), method = "loess", se = T,
              alpha = .45) +
  theme_bw() +
  xlim(0, 5100) +
  theme(strip.background = element_blank(),
        strip.text.y = element_blank(),
        text = element_text(size = 20)) + 
  theme(legend.key.size = unit(1.5, 'cm'),
        legend.text = element_text(size=20))+ 
  xlab("Time (s.)") + ylab("Entropy")+
  scale_color_manual(values = c("#808000", "#FFC200")) +
  scale_fill_manual(values = c("#808000", "#FFC200"))


# Cumulative acceleration
d <- na.omit(bpmDF) %>% 
  group_by(Meeting) %>%
  mutate(acumulada = cumsum(Acceleration)) %>% 
  ggplot(data = ., aes(x=Time, y=acumulada, color = Type, group = Meeting)) + 
#  geom_point(alpha = 0.05)+
  geom_smooth(aes(color = Type, fill = Type), method = "loess", se = T,
              alpha = .45)+
  theme_bw() + xlim(0, 5100)+
  theme(strip.background = element_blank(),
        strip.text.y = element_blank(),
        text = element_text(size = 20)) + 
  theme(legend.key.size = unit(1.5, 'cm'),
        legend.text = element_text(size=20))+ 
  xlab("Time (s.)") + ylab("Cumulative Acceleration")+
  scale_color_manual(values = c("#808000", "#FFC200")) +
  scale_fill_manual(values = c("#808000", "#FFC200"))


h_patch <- (a | b )/(c | d ) + plot_layout(guides = 'collect')
######  Output: Figure 4 ####
h_patch


file_name = paste("./LOESS", ".png", sep = "")
h_patch
ggsave(file_name, h_patch, device = "tiff",
       width=20, height=15, dpi=150, units="in", compression = "lzw")

#


# cumacc <- na.omit(bpmDF) %>% 
#   group_by(Meeting) %>%
#   mutate(acumulada = cumsum(Acceleration)) %>% 
#   ggplot(data = ., aes(x=Time, y=acumulada)) + 
# #  geom_point(color = paleta[1], alpha = 0.40)+ 
#   geom_line(size = 0.5, color = paleta[1])+
#   facet_wrap(~Meeting, nrow = 3)+
#   theme_bw() + xlim(0, 5100)+
#   theme(strip.background = element_blank(),
#         strip.text.y = element_blank(),
#         text = element_text(size = 20)) + 
#   theme(legend.key.size = unit(1.5, 'cm'),
#         legend.text = element_text(size=20))+ 
#   scale_colour_manual(values = paleta[1]) +
#   xlab("Time (s.)") + ylab("Cumulative Acceleration")
# 
# file_name = paste("./.suppAcc", ".png", sep = "")
# cumacc
# ggsave(file_name, cumacc, device = "tiff", 
#        width=20, height=20, dpi=150, units="in", compression = "lzw")
# 

  