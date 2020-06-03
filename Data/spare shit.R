##### Analysis across data sets #####

setwd("~/Documents/Imperial/Morphological Character Hierarchies/Data")
library(ape)
library(ggplot2)
library(ggpubr)
library(cowplot)
library(tidyr)
library(dplyr)
library(RColorBrewer)

### Asher & Hoftreiter, 2006
# paup
paupAH2006 <- read.nexus("Asher & Hofreiter, 2006/Paup-A&H2006.tre")
paupAH2006 <- consensus(paupAH2006)
# morphy
morphyAH2006 <- read.nexus("Asher & Hofreiter, 2006/Morphy-A&H2006.tre")
morphyAH2006 <- consensus(morphyAH2006)
# ana
anaAH2006 <- read.nexus("Asher & Hofreiter, 2006/Anatrees-A&H2006.txt.nex")
anaAH2006 <- consensus(anaAH2006)
## max nodes
maxnodesAH2006 <- length(paupAH2006$tip.label)-1


### Atkinson, 2019
# paup
paupA2019 <- read.nexus("Atkinson, 2019/Paup-A2019.tre")
paupA2019 <- consensus(paupA2019)
# morphy
morphyA2019 <- read.nexus("Atkinson, 2019/Morphy-A2019.tre")
morphyA2019 <- consensus(morphyA2019)
# ana
## max nodes
maxnodesA2019 <- length(paupA2019$tip.label)-1

### CPRS, 2019
# paup
paupCPRS2019 <- read.nexus("Chani-Posse & Ramirez-Salamanca, 2019/Paupcontree-CPRS2019.tre")
# morphy
morphyCPRS2019 <- read.nexus("Chani-Posse & Ramirez-Salamanca, 2019/Morphycontree-CPRS2019.tre")
# ana
## max nodes
maxnodesCPRS2019 <- length(paupCPRS2019$tip.label)-1


### Claeson et al, 2013
# paup
paupC2013 <- read.nexus("Claeson et al, 2013/Paupcontree-C2013.tre")
# morphy
morphyC2013 <- read.nexus("Claeson et al, 2013/Morphycontree-C2013.tre")
# ana
anaC2013 <- read.nexus("Claeson et al, 2013/Anatrees-C2013.nex")
anaC2013 <- consensus(anaC2013)
## max nodes
maxnodesC2013 <- length(paupC2013$tip.label)-1


### Cloutier et al, 2020
# paup
paupC2020 <- read.nexus("Cloutier et al, 2020/Paupcontree-C2020.tre")
# morphy
morphyC2020 <- read.nexus("Cloutier et al, 2020/Morphycontree-C2020.tre")
# ana
anaC2020 <- read.nexus("Cloutier et al, 2020/Anatrees-C2020.nex")
anaC2020 <- consensus(anaC2020)
## max nodes
maxnodesC2020 <- length(paupC2020$tip.label)-1


### Davesne et al, 2016
# paup
paupD2016 <- read.nexus("Davesne et al, 2016/Paup-D2016.tre")
paupD2016 <- consensus(paupD2016)  
# morphy
morphyD2016 <- read.nexus("Davesne et al, 2016/Morphy-D2016.tre")
morphyD2016 <- consensus(morphyD2016)  
# ana
anaD2016 <- read.nexus("Davesne et al, 2016/Anatrees-D2016.nex")
anaD2016 <- consensus(anaD2016)  
## max nodes
maxnodesD2016 <- length(paupD2016$tip.label)-1


### de Lavigerie et al, 2020
# paup
paupdL2020 <- read.nexus("de Lavigerie et al, 2020/Paupcontree-dL2020.tre")
# morphy
morphydL2020 <- read.nexus("de Lavigerie et al, 2020/Morphycontree-dL2020.tre")
# ana
## max nodes
maxnodesdL2020 <- length(paupdL2020$tip.label)-1


### Evans et al, 2014
# paup
paupE2014 <- read.nexus("Evans et al, 2014/Paupcontree-E2014.tre")
# morphy
morphyE2014 <- read.nexus("Evans et al, 2014/Morphycontree-E2014.tre")
# ana
## max nodes
maxnodesE2014 <- length(paupE2014$tip.label)-1


### Fordyce & Marx, 2012
# paup
paupFM2012 <- read.nexus("Fordyce & Marx, 2012/Paupcontree-FM2012.tre")
# morphy
morphyFM2012 <- read.nexus("Fordyce & Marx, 2012/Morphycontree-FM2020.tre")
# ana
## max nodes
maxnodesFM2012 <- length(paupFM2012$tip.label)-1


### Takahashi, 2003
# paup
paupT2003 <- read.nexus("Takahashi, 2003/Paupcontree-T2003.tre")
# morphy
morphyT2003 <- read.nexus("Takahashi, 2003/Morphycontree-T2003.tre")
# ana
## max nodes
maxnodesT2003 <- length(paupT2003$tip.label)-1


### proportion of resolved nodes
# paup
propnodepaupAH2006 <- paupAH2006$Nnode/maxnodesAH2006
propnodepaupA2019 <- paupA2019$Nnode/maxnodesA2019
propnodepaupCPRS2019 <- paupCPRS2019$Nnode/maxnodesCPRS2019
propnodepaupC2013 <- paupC2013$Nnode/maxnodesC2013
propnodepaupC2020 <- paupC2020$Nnode/maxnodesC2020
propnodepaupD2016 <- paupD2016$Nnode/maxnodesD2016
propnodepaupdL2020 <- paupdL2020$Nnode/maxnodesdL2020
propnodepaupE2014 <- paupE2014$Nnode/maxnodesE2014
propnodepaupFM2012 <- paupFM2012$Nnode/maxnodesFM2012
propnodepaupT2003 <- paupT2003$Nnode/maxnodesT2003

# morphy
propnodemorphyAH2006 <- morphyAH2006$Nnode/maxnodesAH2006
propnodemorphyA2019 <- morphyA2019$Nnode/maxnodesA2019
propnodemorphyCPRS2019 <- morphyCPRS2019$Nnode/maxnodesCPRS2019
propnodemorphyC2013 <- morphyC2013$Nnode/maxnodesC2013
propnodemorphyC2020 <- morphyC2020$Nnode/maxnodesC2020
propnodemorphyD2016 <- morphyD2016$Nnode/maxnodesD2016
propnodemorphydL2020 <- morphydL2020$Nnode/maxnodesdL2020
propnodemorphyE2014 <- morphyE2014$Nnode/maxnodesE2014
propnodemorphyFM2012 <- morphyFM2012$Nnode/maxnodesFM2012
propnodemorphyT2003 <- morphyT2003$Nnode/maxnodesT2003

# ana
propnodeanaAH2006 <- anaAH2006$Nnode/maxnodesAH2006
propnodeanaC2013 <- anaC2013$Nnode/maxnodesC2013
propnodeanaC2020 <- anaC2020$Nnode/maxnodesC2020
propnodeanaD2016 <- anaD2016$Nnode/maxnodesD2016

#### Does one programme consistently produce contrees with greater resolutions? ####
## study data frame
study <- c("AH2006","AH2006","AH2006",
           "C2013","C2013","C2013",
           "C2020","C2020","C2020",
           "D2016","D2016","D2016")
df <- data.frame(study)
programme <- c("paup","morphy","anagallis",
               "paup","morphy","anagallis",
               "paup","morphy","anagallis",
               "paup","morphy","anagallis")
df <- cbind(df, programme)


## add resolution to data frame by study
# paup | morphy | ana
resAH2006 <- c(propnodepaupAH2006, propnodemorphyAH2006, propnodeanaAH2006)
resC2013 <- c(propnodepaupC2013, propnodemorphyC2013, propnodeanaC2013)
resC2020 <- c(propnodepaupC2020, propnodemorphyC2020, propnodeanaC2020)
resD2014 <- c(propnodepaupD2016, propnodemorphyD2016, propnodeanaD2016)

resolution <- c(resAH2006,resC2013,resC2020,resD2014)
df <- cbind(df, resolution)

# x-axis study
# y-axis resolution
# colour geom by study

pdf("Figures/Study resolution.pdf", 5.83, 8.27)
ggplot(df, aes(programme, resolution, colour = study, fill = study)) + 
  geom_col(position = "dodge") + theme_pubr()
dev.off()


#### As resolution increases, does Morphy or Anagallis move further away from Paup? ####
# as resolution increases does distance from paup increase?
# x-axis resolution
# y-axis metric
# colour geom by study or pair

df2 <- read.csv("Resolution by distance.csv")

tcres <- ggplot(df2) +
  geom_point(aes(Resolution, Tree.contradiction, colour = Programme)) +
  labs(y = "Tree contradiction from PAUP") +
  theme_pubr()

rfres <- ggplot(df2) +
  geom_point(aes(Resolution, Robinson.Foulds.distance, colour = Programme)) +
  labs(y = "Robinson-Foulds distance from PAUP") +
  theme_pubr()

qdres <- ggplot(df2) +
  geom_point(aes(Resolution, Quartet.distance, colour = Programme)) +
  labs(y = "Quartet distance from PAUP") +
  theme_pubr() + theme(legend.box = "horizontal", legend.position = "top")

distresplots <- plot_grid(tcres + theme(legend.position = "none"),
                          rfres + theme(legend.position = "none"),
                          qdres + theme(legend.position = "none"),
                          labels = "AUTO", nrow = 1, ncol = 3)
distreslegend <- get_legend(qdres +
                              theme(legend.position = "bottom")+
                              guides(colour = guide_legend(nrow = 1)))

pdf("Figures/resolution x distance.pdf", 8, 5.83)
plot_grid(distresplots, distreslegend, ncol = 1, rel_heights = c(1,0.1))
dev.off()

### not a great graph though

#### Does accounting for inapplicability increase resolution? ####
inappcount <- read.csv("Inapp Count.csv")
inappcount$paupmissing <- inappcount$Inapplicable + inappcount$Missing
inappcount$proppaupmissing <- inappcount$paupmissing/inappcount$Total
inappcount$propinapp <- inappcount$Inapplicable/inappcount$Total

paupres <- c(propnodepaupAH2006,
             propnodepaupA2019,
             propnodepaupCPRS2019,
             propnodepaupC2013,
             propnodepaupC2020,
             propnodepaupD2016,
             propnodepaupdL2020,
             propnodepaupE2014,
             propnodepaupFM2012,
             propnodepaupT2003)
inappcount$paupres <- paupres


morphyres <- c(propnodemorphyAH2006,
               propnodemorphyA2019,
               propnodemorphyCPRS2019,
               propnodemorphyC2013,
               propnodemorphyC2020,
               propnodemorphyD2016,
               propnodemorphydL2020,
               propnodemorphyE2014,
               propnodemorphyFM2012,
               propnodemorphyT2003)
inappcount$morphyres <- morphyres


ggplot(inappcount) +
  geom_point(aes(propinapp, morphyres, colour = Study)) +
  labs(x = "Proportion of inapplicable data", y = "Resolution", title = "Morphy") +
  theme_pubr()

ggplot(inappcount) +
  geom_point(aes(proppaupmissing, paupres, colour = Study)) +
  labs(x = "Proportion of missing data", y = "Resolution", title = "PAUP") +
  theme_pubr()


ggplot(inappcount) +
  geom_col(aes(Study, paupres))

ggplot(inappcount) +
  geom_col(aes(Study, morphyres))

barplotdf <- subset(inappcount, select = c("Study", "paupres", "morphyres"))
barplotdf <- pivot_longer(barplotdf, cols = c(paupres, morphyres), names_to = "Programme", values_to = "Resolution")

ggplot(barplotdf2) +
  geom_col(aes(Study, Resolution, fill = Programme), position = "dodge") +
  theme_pubr() +
  fill_palette(palette = "Set2")



