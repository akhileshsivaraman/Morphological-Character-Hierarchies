##### Normalised Distances by Study #####

setwd("~/Documents/Imperial/Morphological Character Hierarchies/Data")
library(ape)
library(ggplot2)
library(ggpubr)
library(cowplot)
library(tidyr)
library(dplyr)
library(RColorBrewer)


# AH2006
AH2006 <- read.csv("Asher & Hofreiter, 2006/normdataAH2006.csv")
AH2006$Study <- c("AH2006","AH2006","AH2006")

# A2019
A2019 <- read.csv("Atkinson, 2019/normdataA2019.csv")
A2019$Study <- c("A2019","A2019","A2019")

# CPRS2019
CPRS2019 <- read.csv("Chani-Posse & Ramirez-Salamanca, 2019/normdataCPRS2019.csv")
CPRS2019$Study <- c("CPRS2019")

# C2013
C2013 <- read.csv("Claeson et al, 2013/normdataC2013.csv")
C2013$Study <- c("C2013","C2013","C2013")

# C2020
C2020 <- read.csv("Cloutier et al, 2020/normdataC2020.csv")
C2020$Study <- c("C2020","C2020","C2020")

# D2016
D2016 <- read.csv("Davesne et al, 2016/normdataD2016.csv")
D2016$Study <- c("D2016","D2016","D2016")

# dL2020
dL2020 <- read.csv("de Lavigerie et al, 2020/normdatadL2020.csv")
dL2020$Study <- c("dL2020","dL2020","dL2020")

# E2014
E2014 <- read.csv("Evans et al, 2014/normdataE2014.csv")
E2014$Study <- c("E2014","E2014","E2014")

# FM2012
FM2012 <- read.csv("Fordyce & Marx, 2012/normdataFM2012.csv")
FM2012$Study <- c("FM2012","FM2012","FM2012")

# L2020
L2020 <- read.csv("Lehtonen, 2020/normdataL2020.csv")
L2020$Study <- c("L2020","L2020","L2020")

# L2015
L2015 <- read.csv("Lihoreau et al, 2015/normdataL2015.csv")
L2015$Study <- c("L2015")

# R2014
R2014 <- read.csv("Rose et al, 2014/normdataR2014.csv")
R2014$Study <- c("R2014")

# T2003
T2003 <- read.csv("Takahashi, 2003/normdataT2003.csv")
T2003$Study <- c("T2003","T2003","T2003")


#### distance by study ####
# x-axis = study
# y-axis = distance

df <- rbind(AH2006, A2019, CPRS2019, C2013, C2020, D2016, dL2020, E2014, FM2012, L2020, L2015, R2014, T2003)

palette <- wesanderson::wes_palette(name = "Chevalier1")
palette <- c(palette, "#25291C", "#E6E49F")
# PAUP* = palette[1]
# Morphy = palette[2]
# Anagallis = palette[3]
# PAUP* vs Morphy = palette[4]
# PAUP* vs Anagallis = palette[5]
# Morphy vs Anagallis = palette[6]


### RF distance
RF <- ggplot(df) +
  geom_pointrange(aes(x = Study, y = Robinson.Foulds.Distance, ymin = Robinson.Foulds.Distance.LQ, 
                      ymax = Robinson.Foulds.Distance.UQ, colour = X), position = position_dodge(width = 0.6), size = 0.6) +
  theme_pubr() +
  color_palette(palette = c(palette[6],palette[5],palette[4])) +
  labs(x = "Study", y = "Normalised Robinson-Foulds distance") +
  guides(colour = guide_legend(title = "Pairwise comparison")) +
  theme(legend.text = element_text(size = 12), axis.text.x = element_text(size=9.5))


### Quartet distance
Q <- ggplot(df) +
  geom_pointrange(aes(x = Study, y = Quartet.Distances, ymin = Quartet.Distances.LQ, 
                      ymax = Quartet.Distances.UQ, colour = X), position = position_dodge(width = 0.6), size = 0.6) +
  theme_classic2() +
  color_palette(palette = c(palette[6],palette[5],palette[4])) +
  labs(x = "Study", y = "Normalised quartet distance") +
  guides(colour = guide_legend(title = "Pairwise comparison")) +
  theme(legend.text = element_text(size = 12), axis.text.x = element_text(size=9.5), 
        legend.box = "horizontal", legend.position = "top")


plots <- plot_grid(RF + theme(legend.position = "none"),
                   Q + theme(legend.position = "none"),
                   labels = "AUTO", nrow = 2, ncol = 1)

legend <- get_legend(Q +
                       theme(legend.position = "top") +
                       guides(colour = guide_legend(ncol = 3, title = "Pairwise comparison")))

pdf("Figures/distances by study.pdf", 8.27, 11.5)
plot_grid(legend, plots, ncol = 1, rel_heights = c(0.1,1))
dev.off()


write.csv(df, "Normalised distances for all studies.csv")









