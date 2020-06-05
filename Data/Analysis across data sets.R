##### Analysis across data sets #####

setwd("~/Documents/Imperial/Morphological Character Hierarchies/Data")
library(ape)
library(ggplot2)
library(ggpubr)
library(cowplot)
library(tidyr)
library(dplyr)
library(data.table)
library(RColorBrewer)

palette <- wesanderson::wes_palette(name = "Chevalier1")
palette <- c(palette, "#25291C", "#E6E49F")
# PAUP* = palette[1]
# Morphy = palette[2]
# Anagallis = palette[3]
# PAUP* vs Morphy = palette[4]
# PAUP* vs Anagallis = palette[5]
# Morphy vs Anagallis = palette[6]

### Asher & Hoftreiter, 2006
# paup
paupAH2006 <- consensus(read.nexus("Asher & Hofreiter, 2006/Paup-A&H2006.tre"))
# morphy
morphyAH2006 <- consensus(read.nexus("Asher & Hofreiter, 2006/Morphy-A&H2006.tre"))
# ana
anaAH2006 <- consensus(read.nexus("Asher & Hofreiter, 2006/Anatrees-A&H2006.txt.nex"))


### Atkinson, 2019
# paup
paupA2019 <- read.nexus("Atkinson, 2019/Paupcontree-A2019.tre")
# morphy
morphyA2019 <- consensus(read.nexus("Atkinson, 2019/Morphytrees-A2019.tre"))
# ana
anaA2019 <- consensus(read.nexus("Atkinson, 2019/Anatrees-A2019.nex"))

### CPRS, 2019
# paup
paupCPRS2019 <- read.nexus("Chani-Posse & Ramirez-Salamanca, 2019/Paupcontree-CPRS2019.tre")
# morphy
morphyCPRS2019 <- read.nexus("Chani-Posse & Ramirez-Salamanca, 2019/Morphycontree-CPRS2019.tre")


### Claeson et al, 2013
# paup
paupC2013 <- read.nexus("Claeson et al, 2013/Paupcontree-C2013.tre")
# morphy
morphyC2013 <- read.nexus("Claeson et al, 2013/Morphycontree-C2013.tre")
# ana
anaC2013 <- consensus(read.nexus("Claeson et al, 2013/Anatrees-C2013.nex"))


### Cloutier et al, 2020
# paup
paupC2020 <- read.nexus("Cloutier et al, 2020/Paupcontree-C2020.tre")
# morphy
morphyC2020 <- read.nexus("Cloutier et al, 2020/Morphycontree-C2020.tre")
# ana
anaC2020 <- consensus(read.nexus("Cloutier et al, 2020/Anatrees-C2020.nex"))


### Davesne et al, 2016
# paup
paupD2016 <- consensus(read.nexus("Davesne et al, 2016/Paup-D2016.tre"))
# morphy
morphyD2016 <- consensus(read.nexus("Davesne et al, 2016/Morphy-D2016.tre"))
# ana
anaD2016 <- consensus(read.nexus("Davesne et al, 2016/Anatrees-D2016.nex"))


### de Lavigerie et al, 2020
# paup
paupdL2020 <- read.nexus("de Lavigerie et al, 2020/Paupcontree-dL2020.tre")
# morphy
morphydL2020 <- read.nexus("de Lavigerie et al, 2020/Morphycontree-dL2020.tre")
# ana
anadL2020 <- consensus(read.nexus("de Lavigerie et al, 2020/Anatrees-dL2020.txt.nex"))

### Evans et al, 2014
# paup
paupE2014 <- read.nexus("Evans et al, 2014/Paupcontree-E2014.tre")
# morphy
morphyE2014 <- read.nexus("Evans et al, 2014/Morphycontree-E2014.tre")
# ana
anaE2014 <- consensus(read.nexus("Evans et al, 2014/Anatrees-E2014.txt.nex"))

### Fordyce & Marx, 2012
# paup
paupFM2012 <- read.nexus("Fordyce & Marx, 2012/Paupcontree-FM2012.tre")
# morphy
morphyFM2012 <- read.nexus("Fordyce & Marx, 2012/Morphycontree-FM2020.tre")
# ana
anaFM2012 <- consensus(read.nexus("Fordyce & Marx, 2012/Anatrees-FM2012.txt.nex"))

### Lehtonen et al, 2020
# paup
paupL2020 <- read.nexus("Lehtonen, 2020/Paupcontree-L2020.tre")
# morphy
morphyL2020 <- read.nexus("Lehtonen, 2020/Morphycontree-L2020.tre")
# ana
anaL2020 <- consensus(read.nexus("Lehtonen, 2020/Anatrees-L2020.nex"))

### Lihoreau et al, 2015
# paup
paupL2015 <- read.nexus("Lihoreau et al, 2015/Paupcontree-L2015.tre")
# morphy
morphyL2015 <- read.nexus("Lihoreau et al, 2015/Morphycontree-L2015.tre")

### Rose, 2014
# paup
paupR2014 <- read.nexus("Rose et al, 2014/Paupcontree-R2014.tre")
# morphy
morphyR2014 <- read.nexus("Rose et al, 2014/Morphycontree-R2014.tre")


### Takahashi, 2003
# paup
paupT2003 <- read.nexus("Takahashi, 2003/Paupcontree-T2003.tre")
# morphy
morphyT2003 <- read.nexus("Takahashi, 2003/Morphycontree-T2003.tre")
# ana
anaT2003 <- consensus(read.nexus("Takahashi, 2003/Anatrees-T2003.nex"))


#### Resolution metric ####
### max nodes
maxnodesAH2006 <- length(paupAH2006$tip.label)-1
maxnodesA2019 <- length(paupA2019$tip.label)-1
maxnodesCPRS2019 <- length(paupCPRS2019$tip.label)-1
maxnodesC2013 <- length(paupC2013$tip.label)-1
maxnodesC2020 <- length(paupC2020$tip.label)-1
maxnodesD2016 <- length(paupD2016$tip.label)-1
maxnodesdL2020 <- length(paupdL2020$tip.label)-1
maxnodesE2014 <- length(paupE2014$tip.label)-1
maxnodesFM2012 <- length(paupFM2012$tip.label)-1
maxnodesL2020 <- length(paupL2020$tip.label)-1
maxnodesL2015 <- length(paupL2015$tip.label)-1
maxnodesR2014 <- length(paupR2014$tip.label)-1
maxnodesT2003 <- length(paupT2003$tip.label)-1


### propotion of resolved nodes = resolution metric
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
propnodepaupL2020 <- paupL2020$Nnode/maxnodesL2020
propnodepaupL2015 <- paupL2015$Nnode/maxnodesL2015
propnodepaupR2014 <- paupR2014$Nnode/maxnodesR2014
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
propnodemorphyL2020 <- morphyL2020$Nnode/maxnodesL2020
propnodemorphyL2015 <- morphyL2015$Nnode/maxnodesL2015
propnodemorphyR2014 <- morphyR2014$Nnode/maxnodesR2014
propnodemorphyT2003 <- morphyT2003$Nnode/maxnodesT2003


# ana
propnodeanaAH2006 <- anaAH2006$Nnode/maxnodesAH2006
propnodeanaA2019 <- anaA2019$Nnode/maxnodesA2019
propnodeanaC2013 <- anaC2013$Nnode/maxnodesC2013
propnodeanaC2020 <- anaC2020$Nnode/maxnodesC2020
propnodeanaD2016 <- anaD2016$Nnode/maxnodesD2016
propnodeanadL2020 <- anadL2020$Nnode/maxnodesdL2020
propnodeanaE2014 <- anaE2014$Nnode/maxnodesE2014
propnodeanaFM2012 <- anaFM2012$Nnode/maxnodesFM2012
propnodeanaL2020 <- anaL2020$Nnode/maxnodesL2020
propnodeanaT2003 <- anaT2003$Nnode/maxnodesT2003



##### is there a correlation between inapplicable count/proportion and the magnitude of the distance from PAUP? ####
inappdist <- read.csv("Resolution by distance.csv")

shapiro.test(inappdist$Robinson.Foulds.distance) # normal
shapiro.test(inappdist$Quartet.distance) # normal

### inapplicable proportion
shapiro.test(inappdist$Proportion.of.inapplicable) # non-normal
# use non-parametric test
corInappPropXRF <- lapply(split(inappdist, inappdist$Programme), function(x)
  cor.test(x$Proportion.of.inapplicable, x$Robinson.Foulds.distance, method = "spearman"))
corInappPropXRF
# Anagallis: rho = -0.05454545, S = 174, p-value = 0.8916
# Morphy: rho = -0.005494505, S = 366, p-value = 0.9928

corInappPropXQu <- lapply(split(inappdist, inappdist$Programme), function(x)
  cor.test(x$Proportion.of.inapplicable, x$Quartet.distance, method = "spearman"))
corInappPropXQu
# Anagallis: rho = -0.006060606, S = 356, p-value = 0.9494
# Morphy: rho = 0.02197802, S = 160, p-value = 0.9457


plotInappPropXRF <- ggplot(inappdist) +
  geom_point(aes(Proportion.of.inapplicable, Robinson.Foulds.distance, colour = Programme), size = 2.5) +
  labs(x = "Proportion of 
inapplicable data", y = "Normalised Robinson-Foulds distance") +
  theme_pubr() +
  color_palette(palette = c(palette[3], palette[2]))

plotInappPropXQu <- ggplot(inappdist) +
  geom_point(aes(Proportion.of.inapplicable, Quartet.distance, colour = Programme), size = 2.5) +
  labs(x = "Proportion of 
inapplicable data", y = "Normalised quartet distance") +
  theme_pubr() +
  color_palette(palette = c(palette[3], palette[2])) +
  theme(legend.box = "horizontal", legend.position = "top")


plotsInappPropxDistance <- plot_grid(plotInappPropXRF + theme(legend.position = "none", axis.title.x = element_text(size = 11), axis.text = element_text(size = 10)),
                            plotInappPropXQu + theme(legend.position = "none", axis.title.x = element_text(size = 11), axis.text = element_text(size = 10)),
                            labels = "AUTO", nrow = 1, ncol = 2)
legendInappPropxDistance <- get_legend(plotInappPropXQu +
                                theme(legend.position = "bottom")+
                                guides(colour = guide_legend(nrow = 1)))

pdf("Figures/is there a correlation between inapp proportion and magnitude of the distance from paup.pdf", 8, 5.83)
plot_grid(plotsInappPropxDistance, legendInappPropxDistance, ncol = 1, rel_heights = c(1,0.1))
dev.off()


### inapplicable count
shapiro.test(inappdist$Inapplicable.count) # non-normal
# use non-parametric test
corInappCountXRF <- lapply(split(inappdist, inappdist$Programme), function(x)
  cor.test(x$Inapplicable.count, x$Robinson.Foulds.distance, method = "spearman"))
corInappCountXRF
# Anagallis: rho = -0.1878788, S = 196, p-value = 0.6076
# Morphy: rho = -0.5, S = 546, p-value = 0.08498

corInappCountXQu <- lapply(split(inappdist, inappdist$Programme), function(x)
  cor.test(x$Inapplicable.count, x$Quartet.distance, method = "spearman"))
corInappCountXQu
# Anagallis: rho = -0.2484848, S = 206, p-value = 0.4916
# Morphy: rho = -0.4835165, S = 540, p-value = 0.09704


plotInappCountXRF <- ggplot(inappdist) +
  geom_point(aes(Inapplicable.count, Robinson.Foulds.distance, colour = Programme), size = 2.5) +
  labs(x = "Inapplicable count", y = "Normalised Robinson-Foulds distance") +
  theme_pubr() +
  color_palette(palette = c(palette[3], palette[2]))

plotInappCountXQu <- ggplot(inappdist) +
  geom_point(aes(Inapplicable.count, Quartet.distance, colour = Programme), size = 2.5) +
  labs(x = "Inapplicable count", y = "Normalised quartet distance") +
  theme_pubr() +
  color_palette(palette = c(palette[3], palette[2])) +
  theme(legend.box = "horizontal", legend.position = "top")


plotsInappCountxDistance <- plot_grid(plotInappCountXRF + theme(legend.position = "none", axis.title.x = element_text(size = 11), axis.text = element_text(size = 10)),
                                     plotInappCountXQu + theme(legend.position = "none", axis.title.x = element_text(size = 11), axis.text = element_text(size = 10)),
                                     labels = "AUTO", nrow = 1, ncol = 2)
legendInappCountxDistance <- get_legend(plotInappCountXQu +
                                         theme(legend.position = "bottom")+
                                         guides(colour = guide_legend(nrow = 1)))

pdf("Figures/is there a correlation between inapp count and magnitude of the distance from paup.pdf", 8, 5.83)
plot_grid(plotsInappCountxDistance, legendInappCountxDistance, ncol = 1, rel_heights = c(1,0.1))
dev.off()

allplots <- plot_grid(plotInappPropXRF + theme(legend.position = "none", axis.title.x = element_text(size = 11), axis.text = element_text(size = 10)),
          plotInappPropXQu + theme(legend.position = "none", axis.title.x = element_text(size = 11), axis.text = element_text(size = 10)),
          plotInappCountXRF + theme(legend.position = "none", axis.title.x = element_text(size = 11), axis.text = element_text(size = 10)),
          plotInappCountXQu + theme(legend.position = "none", axis.title.x = element_text(size = 11), axis.text = element_text(size = 10)),
          labels = "AUTO", nrow = 2, ncol = 2)
plot_grid(allplots, legendInappCountxDistance, ncol = 1, rel_heights = c(1,0.1))


##### is there a difference betweeen inferring inapplicability and explicitly stating inapplicability? ####
### which one is closer to PAUP?
# normalised distance by study graph




#### is there a correlation between resolution and distance from paup? ####
# x-axis = resolution
# y-axis = distance from paup
resbydist <- read.csv("Resolution by distance.csv")

## correlation test
# is resolution normally distributed?
tapply(resbydist$Resolution, resbydist$Programme, shapiro.test) # yes

# is distance from paup normally distributed?
tapply(resbydist$Robinson.Foulds.distance, resbydist$Programme, shapiro.test) # yes
tapply(resbydist$Quartet.distance, resbydist$Programme, shapiro.test) # yes

# correlation tests
corresxrf <- lapply(split(resbydist, resbydist$Programme), function(x)
  cor.test(x$Resolution, x$Robinson.Foulds.distance))
corresxrf
# Anagallis: r = -0.6523977, t = -2.4348, df = 8, p-value = 0.0409
# Morphy: r = -0.637404, t = -2.7436, df = 11, p-value = 0.01911

corresxqd <- lapply(split(resbydist, resbydist$Programme), function(x)
  cor.test(x$Resolution, x$Quartet.distance, ))
corresxqd
# Anagallis: r = -0.671958, t = -2.5663, df = 8, p-value = 0.03332
# Morphy: r = -0.3025855, t = -1.0529, df = 11, p-value = 0.315


rfres <- ggplot(resbydist) +
  geom_point(aes(Resolution, Robinson.Foulds.distance, colour = Programme), size = 2) +
  geom_smooth(aes(Resolution, Robinson.Foulds.distance, colour = Programme), method = "lm",
              fullrange=T, se = F, size = 0.75) +
  labs(y = "Normalised Robinson-Foulds distance from PAUP*") +
  theme_pubr() +
  color_palette(palette = c(palette[3], palette[2]))

qdres <- ggplot(resbydist) +
  geom_point(aes(Resolution, Quartet.distance, colour = Programme), size = 2) +
  geom_smooth(aes(Resolution, Quartet.distance, colour = Programme), method = "lm",
              fullrange=T, se = F, size = 0.75) +
  labs(y = "Normlised quartet distance from PAUP*") +
  theme_pubr() + 
  color_palette(palette = c(palette[3], palette[2])) +
  theme(legend.box = "horizontal", legend.position = "top")

distresplots <- plot_grid(rfres + theme(legend.position = "none"),
                          qdres + theme(legend.position = "none"),
                          labels = "AUTO", nrow = 1, ncol = 2)
distreslegend <- get_legend(qdres +
                              theme(legend.position = "bottom")+
                              guides(colour = guide_legend(nrow = 1)))

pdf("Figures/is there a correlation between resolution and distance from paup.pdf", 8, 5.83)
plot_grid(distresplots, distreslegend, ncol = 1, rel_heights = c(1,0.1))
dev.off()


#### is there a correlation between inapp proportion/count and resolution? ####
inappres <- read.csv("Resolution by distance.csv")


## correlation tests
# is resolution normally distributed?
tapply(inappres$Resolution, inappres$Programme, shapiro.test) # yes

# are inapp data normally distributed?
tapply(inappres$Inapplicable.count, inappres$Programme, shapiro.test) # no
tapply(inappres$Proportion.of.inapplicable, inappres$Programme, shapiro.test) # no

corinappcxres <- lapply(split(resbydist, resbydist$Programme), function(x)
  cor.test(x$Resolution, x$Inapplicable.count))
corinappcxres
# Anagallis: r = 0.463201, t = 1.4783, df = 8, p-value = 0.1776
# Morphy: r = 0.3016198, t = 1.0492, df = 11, p-value = 0.3166

corinapppropxres <- lapply(split(resbydist, resbydist$Programme), function(x)
  cor.test(x$Resolution, x$Proportion.of.inapplicable))
corinapppropxres
# Anagallis: r = -0.3134446, t = -0.9336, df = 8, p-value = 0.3778
# Morphy: r = -0.1606066, t = -0.53968, df = 11, p-value = 0.6002

inappcres <- ggplot(inappres) +
  geom_point(aes(Inapplicable.count, Resolution, colour = Programme), size = 2) +
  geom_smooth(aes(Inapplicable.count, Resolution, colour = Programme), method = "lm",
              fullrange=T, se = F, size = 0.75) +
  labs(x = "Inapplicable count") +
  theme_pubr() +
  color_palette(palette = c(palette[3], palette[2]))

inapppropres <- ggplot(inappres) +
  geom_point(aes(Proportion.of.inapplicable, Resolution, colour = Programme), size = 2) +
  geom_smooth(aes(Proportion.of.inapplicable, Resolution, colour = Programme), method = "lm",
              fullrange=T, se = F, size = 0.75) +
  labs(x = "Proportion of inapplicable data") +
  theme_pubr() +
  color_palette(palette = c(palette[3], palette[2]))

inappresplots <- plot_grid(inappcres + theme(legend.position = "none"),
                           inapppropres + theme(legend.position = "none"),
                           labels = "AUTO", nrow = 1, ncol = 2)
inappreslegend <- get_legend(inapppropres +
                               theme(legend.position = "bottom")+
                               guides(colour = guide_legend(nrow = 1)))

pdf("Figures/is there a correlation between resolution and inapp.pdf", 8, 5.83)
plot_grid(inappresplots, inappreslegend, ncol = 1, rel_heights = c(1,0.1))
dev.off()





#### is there a difference in the proportion of resolved nodes ####
## barplot
# y-axis = resolution
# x-axis = study, each study with multiple bars, one for each programme

# Study | Resolution | Programme
study <- c("AH2006",
           "A2019",
           "CPRS2019",
           "C2013",
           "C2020",
           "D2016",
           "dL2020",
           "E2014",
           "FM2012",
           "L2020",
           "L2015",
           "R2014",
           "T2003")
paupres <- c(propnodepaupAH2006,
             propnodepaupA2019,
             propnodepaupCPRS2019,
             propnodepaupC2013,
             propnodepaupC2020,
             propnodepaupD2016,
             propnodepaupdL2020,
             propnodepaupE2014,
             propnodepaupFM2012,
             propnodepaupL2020,
             propnodepaupL2015,
             propnodepaupR2014,
             propnodepaupT2003)
morphyres <- c(propnodemorphyAH2006,
               propnodemorphyA2019,
               propnodemorphyCPRS2019,
               propnodemorphyC2013,
               propnodemorphyC2020,
               propnodemorphyD2016,
               propnodemorphydL2020,
               propnodemorphyE2014,
               propnodemorphyFM2012,
               propnodemorphyL2020,
               propnodemorphyL2015,
               propnodemorphyR2014,
               propnodemorphyT2003)

propnodeanaCPRS2019 <- 0
propnodeanaL2015 <- 0
propnodeanaR2014 <- 0

anares <- c(propnodeanaAH2006,
            propnodeanaA2019,
            propnodeanaCPRS2019,
            propnodeanaC2013,
            propnodeanaC2020,
            propnodeanaD2016,
            propnodeanadL2020,
            propnodeanaE2014,
            propnodeanaFM2012,
            propnodeanaL2020,
            propnodeanaL2015,
            propnodeanaR2014,
            propnodeanaT2003)

resbyprogramme <- as.data.frame(study)
resbyprogramme$paupres <- paupres
resbyprogramme$morphyres <- morphyres
resbyprogramme$anares <- anares
colnames(resbyprogramme) <- c("Study", "PAUP*", "Morphy", "Anagallis")
resbyprogramme <- pivot_longer(resbyprogramme, cols = c("PAUP*", "Morphy", "Anagallis"),
                               names_to = "Programme", values_to = "Resolution")


pdf("Figures/is there a difference in the proprotion of resolved nodes.pdf", 11, 8)
ggplot(resbyprogramme) +
  geom_col(aes(Study, Resolution, fill = Programme), position = "dodge") +
  theme_pubr() +
  fill_palette(palette = c(palette[3], palette[2], palette[1]))
dev.off()

## test if they are significantly different? how?



########## extra shit ##########

#### morphy and paup: does resolution vary between programmes? ####
pmresbyprogramme <- read.csv("Inapp Count.csv")
pmresbyprogramme <- subset(pmresbyprogramme, select = c("Study", "Inapplicable"))
allmorphyres <- c(propnodemorphyAH2006,
                  propnodemorphyC2013,
                  propnodemorphyC2020,
                  propnodemorphyD2016,
                  propnodemorphyA2019,
                  propnodemorphyCPRS2019,
                  propnodemorphydL2020,
                  propnodemorphyE2014,
                  propnodemorphyFM2012,
                  propnodemorphyT2003,
                  propnodemorphyL2015,
                  propnodemorphyR2014)
allpaupres <- c(propnodepaupAH2006,
                propnodepaupC2013,
                propnodepaupC2020,
                propnodepaupD2016,
                propnodepaupA2019,
                propnodepaupCPRS2019,
                propnodepaupdL2020,
                propnodepaupE2014,
                propnodepaupFM2012,
                propnodepaupT2003,
                propnodepaupL2015,
                propnodepaupR2014)
pmresbyprogramme <- cbind(pmresbyprogramme, allmorphyres, allpaupres)
colnames(pmresbyprogramme) <- c("Study", "Inapplicable", "Morphy", "PAUP*")

pmresbyprogramme <- pivot_longer(pmresbyprogramme, cols = c("Morphy", "PAUP*"),
                               names_to = "Programme", values_to = "Resolution")

pdf("Figures/does resolution vary between paup and morphy.pdf", 8, 5.83)
ggplot(pmresbyprogramme) +
  geom_col(aes(Study, Resolution, fill = Programme), position = "dodge") +
  theme_pubr() +
  fill_palette(palette = c(palette[2], palette[1]))
dev.off()


#### all: is there a correlation between inapp proportion/count and distance from paup? ####
inappdist <- read.csv("Resolution by distance.csv")

### proportion of inapp data
## correlation test
# distances from paup is distributed normally - for tests check "is there a correlation between resolution and distance from paup"
# proportion of inapp data is not normally distributed - for test check "is there a correlation between inapp and resolution"
# use non-parametric test
corinappxtc <- lapply(split(inappdist, inappdist$Programme), function(x)
  cor.test(x$Proportion.of.inapplicable, x$Tree.contradiction, method = "spearman"))
# Anagallis: rho=-0.4, p=0.75 => non-significant -ve correlation
# Morphy: rho=-0.2, p=0.9167 => non-significant -ve correlation

corinappxrf <- lapply(split(inappdist, inappdist$Programme), function(x)
  cor.test(x$Proportion.of.inapplicable, x$Robinson.Foulds.distance, method = "spearman"))
# Anagallis: rho=-0.4, p=0.75 => non-significant -ve correlation
# Morphy: rho=-0.2, p=0.9167 => non-significant -ve correlation

corinappxqd <- lapply(split(inappdist, inappdist$Programme), function(x)
  cor.test(x$Proportion.of.inapplicable, x$Quartet.distance, method = "spearman"))
# Anagallis: rho=-0.4, p=0.75 => non-significant -ve correlation
# Morphy: rho=-0.2, p=0.9167 => non-significant -ve correlation


# x=prop inapp
# y=distance
inapptc <- ggplot(inappdist) +
  geom_point(aes(Proportion.of.inapplicable, Tree.contradiction, colour = Programme), size = 2) +
  labs(x = "Proportion of inapplicable data", y = "Normalised tree contradiction from PAUP*") +
  scale_x_continuous(limits = c(0,0.05)) +
  theme_pubr() +
  color_palette(palette = c(palette[3], palette[2]))

inapprf <- ggplot(inappdist) +
  geom_point(aes(Proportion.of.inapplicable, Robinson.Foulds.distance, colour = Programme), size = 2) +
  labs(x = "Proportion of inapplicable data", y = "Normalised Robinson-Foulds distance from PAUP*") +
  scale_x_continuous(limits = c(0,0.05)) +
  theme_pubr() +
  color_palette(palette = c(palette[3], palette[2]))

inappqd <- ggplot(inappdist) +
  geom_point(aes(Proportion.of.inapplicable, Quartet.distance, colour = Programme), size = 2) +
  labs(x = "Proportion of inapplicable data", y = "Normalised quartet distance from PAUP*") +
  scale_x_continuous(limits = c(0,0.05)) + 
  theme_pubr() +
  color_palette(palette = c(palette[3], palette[2])) + 
  theme(legend.box = "horizontal", legend.position = "top")

inappdistplots <- plot_grid(inapptc + theme(legend.position = "none", axis.title.x = element_text(size = 11), axis.text = element_text(size = 10)),
                          inapprf + theme(legend.position = "none", axis.title.x = element_text(size = 11), axis.text = element_text(size = 10)),
                          inappqd + theme(legend.position = "none", axis.title.x = element_text(size = 11), axis.text = element_text(size = 10)),
                          labels = "AUTO", nrow = 1, ncol = 3)
inappdistlegend <- get_legend(inappqd +
                              theme(legend.position = "bottom")+
                              guides(colour = guide_legend(nrow = 1)))

pdf("Figures/is there a correlation between inapp proportion and distance from paup?.pdf", 8, 5.83)
plot_grid(inappdistplots, inappdistlegend, ncol = 1, rel_heights = c(1,0.1))
dev.off()


### Inapp count
## correlation test
# is inapp count normally distributed?
shapiro.test(inappdist$Inapplicable.count) # yes
# distances from paup are also normally distributed => Pearson's

corinappcountxtc <- lapply(split(inappdist, inappdist$Programme), function(x)
  cor.test(x$Inapplicable.count, x$Tree.contradiction, method = "pearson"))
# Anagallis: r=-0.8109681, p-value=0.189 => non-significant -ve correlation
# Morphy: r=-0.8578704, p-value=0.1421 => non-significant -ve correlation

corinappcountxrf <- lapply(split(inappdist, inappdist$Programme), function(x)
  cor.test(x$Inapplicable.count, x$Robinson.Foulds.distance, method = "pearson"))
# Anagallis: r=-0.2871362, p-value=0.7129 => non-significant -ve correlation
# Morphy: r=-0.8089681, p-value=0.191 => non-significant -ve correlation

corinappcountxqd <- lapply(split(inappdist, inappdist$Programme), function(x)
  cor.test(x$Inapplicable.count, x$Quartet.distance, method = "pearson"))
# Anagallis: r=-0.5223941, p-value=0.4776 => non-significant -ve correlation
# Morphy: r=-0.9035562, p-value=0.09644 => non-significant -ve correlation


# x=prop count
# y=distance
inappcounttc <- ggplot(inappdist) +
  geom_point(aes(Inapplicable.count, Tree.contradiction, colour = Programme), size = 2) +
  labs(x = "Inapplicable count", y = "Normalised tree contradiction from PAUP*") +
  theme_pubr() +
  color_palette(palette = c(palette[3], palette[2]))

inappcountrf <- ggplot(inappdist) +
  geom_point(aes(Inapplicable.count, Robinson.Foulds.distance, colour = Programme), size = 2) +
  labs(x = "Inapplicable count", y = "Normalised Robinson-Foulds distance from PAUP*") +
  theme_pubr() +
  color_palette(palette = c(palette[3], palette[2]))

inappcountqd <- ggplot(inappdist) +
  geom_point(aes(Inapplicable.count, Quartet.distance, colour = Programme), size = 2) +
  labs(x = "Inapplicable count", y = "Normalised quartet distance from PAUP*") +
  theme_pubr() +
  color_palette(palette = c(palette[3], palette[2])) + 
  theme(legend.box = "horizontal", legend.position = "top")

inappcountdistplots <- plot_grid(inappcounttc + theme(legend.position = "none", axis.title.x = element_text(size = 11), axis.text = element_text(size = 10)),
                            inappcountrf + theme(legend.position = "none", axis.title.x = element_text(size = 11), axis.text = element_text(size = 10)),
                            inappcountqd + theme(legend.position = "none", axis.title.x = element_text(size = 11), axis.text = element_text(size = 10)),
                            labels = "AUTO", nrow = 1, ncol = 3)
inappcountdistlegend <- get_legend(inappcountqd +
                                theme(legend.position = "bottom")+
                                guides(colour = guide_legend(nrow = 1)))

pdf("Figures/is there a correlation between inapp count and distance from paup?.pdf", 8, 5.83)
plot_grid(inappcountdistplots, inappcountdistlegend, ncol = 1, rel_heights = c(1,0.1))
dev.off()



#### morphy and paup: is there a correlation between inapp proportion/count and distance from paup? ####
pminappdist <- read.csv("Inapp Count.csv")
normAH2006 <- read.csv("Asher & Hofreiter, 2006/normdataAH2006.csv")
normAH2006 <- normAH2006[1,]
normA2019 <- read.csv("Atkinson, 2019/normdataA2019.csv")
normCPRS2019 <- read.csv("Chani-Posse & Ramirez-Salamanca, 2019/normdataCPRS2019.csv")
normC2013 <- read.csv("Claeson et al, 2013/normdataC2013.csv")
normC2013 <- normC2013[1,]
normC2020 <- read.csv("Cloutier et al, 2020/normdataC2020.csv")
normC2020 <- normC2020[1,]
normD2016 <- read.csv("Davesne et al, 2016/normdataD2016.csv")
normD2016 <- normD2016[1,]
normdL2020 <- read.csv("de Lavigerie et al, 2020/normdatadL2020.csv")
normE2014 <- read.csv("Evans et al, 2014/normdataE2014.csv")
normFM2012 <- read.csv("Fordyce & Marx, 2012/normdataFM2012.csv")
normT2003 <- read.csv("Takahashi, 2003/normdataT2003.csv")
normL2015 <- read.csv("Lihoreau et al, 2015/normdataL2015.csv")
normR2014 <- read.csv("Rose et al, 2014/normdataR2014.csv")

a <- rbind(normAH2006,
           normA2019,
           normCPRS2019,
           normC2013,
           normC2020,
           normD2016,
           normdL2020,
           normE2014,
           normFM2012,
           normT2003,
           normL2015,
           normR2014)
pminappdist <- cbind(pminappdist, a)


### proportion by distance
## correlation test
cor.test(pminappdist$Proportion.of.inapplicable, pminappdist$Tree.Contradiction, method = "spearman")
# rho=-0.4055944, p=0.1926
cor.test(pminappdist$Proportion.of.inapplicable, pminappdist$Robinson.Foulds.Distance, method = "spearman")
# rho=-0.3846154, p=0.2184
cor.test(pminappdist$Proportion.of.inapplicable, pminappdist$Quartet.Distances, method = "spearman")
# rho=-0.1328671, p=0.6834

## plots
# tree contradiction
propinapptcplot <- ggplot(pminappdist) +
  geom_point(aes(Proportion.of.inapplicable, Tree.Contradiction), colour = palette[2], size = 2) +
  theme_pubr() +
  labs(x = "Proportion of inapplicable data", y = "Normalised tree contradiction from PAUP*")

# RF distance
propinapprfplot <- ggplot(pminappdist) +
  geom_point(aes(Proportion.of.inapplicable, Robinson.Foulds.Distance), colour = palette[2], size = 2) +
  theme_pubr() +
  labs(x = "Proportion of inapplicable data", y = "Normalised Robinson-Foulds distance from PAUP*")

# quartet distance
propinappqdplot <- ggplot(pminappdist) +
  geom_point(aes(Proportion.of.inapplicable, Quartet.Distances), colour = palette[2], size = 2) +
  theme_pubr() +
  labs(x = "Proportion of inapplicable data", y = "Normalised quartet distance from PAUP*")

pdf("Figures/is there a correlation between inapp prop and distance from paup? morphy.pdf", 8, 5.83)
plot_grid(propinapptcplot + theme(legend.position = "none", axis.title.x = element_text(size = 11), axis.text = element_text(size = 10)), 
          propinapprfplot + theme(legend.position = "none", axis.title.x = element_text(size = 11), axis.text = element_text(size = 10)), 
          propinappqdplot + theme(legend.position = "none", axis.title.x = element_text(size = 11), axis.text = element_text(size = 10)), 
          labels = "AUTO", nrow = 1, ncol = 3)
dev.off()


### count by distance
## correlation test
cor.test(pminappdist$Inapplicable, pminappdist$Tree.Contradiction, method = "spearman")
# rho=-0.6013986, p=0.04281 => sig negative correlation
cor.test(pminappdist$Inapplicable, pminappdist$Robinson.Foulds.Distance, method = "spearman")
# rho=-0.5944056, p=0.04575 => sig negative correlation
cor.test(pminappdist$Inapplicable, pminappdist$Quartet.Distances, method = "spearman")
# rho=-0.3146853, p=0.3195

# tree contradiction
countinapptcplot <- ggplot(pminappdist) +
  geom_point(aes(Inapplicable, Tree.Contradiction), colour = palette[2], size = 2) +
  theme_pubr() +
  labs(x = "Inapplicable data count", y = "Normalised tree contradiction from PAUP*")

# RF distance
countinapprfplot <- ggplot(pminappdist) +
  geom_point(aes(Inapplicable, Robinson.Foulds.Distance), colour = palette[2], size = 2) +
  theme_pubr() +
  labs(x = "Inapplicable data count", y = "Normalised Robinson-Foulds distance from PAUP*")

# quartet distance
countinappqdplot <- ggplot(pminappdist) +
  geom_point(aes(Inapplicable, Quartet.Distances), colour = palette[2], size = 2) +
  theme_pubr() +
  labs(x = "Inapplicable data count", y = "Normalised quartet distance from PAUP*")

pdf("Figures/is there a relationship between inapp count and distance from paup? morphy.pdf", 8, 5.83)
plot_grid(countinapptcplot + theme(legend.position = "none", axis.title.x = element_text(size = 11), axis.text = element_text(size = 10)), 
          countinapprfplot + theme(legend.position = "none", axis.title.x = element_text(size = 11), axis.text = element_text(size = 10)), 
          countinappqdplot + theme(legend.position = "none", axis.title.x = element_text(size = 11), axis.text = element_text(size = 10)), 
          labels = "AUTO", nrow = 1, ncol = 3)
dev.off()




#### morphy: is there a correlation between inapp proportion/count and resolution? ####
inappres <- read.csv("Inapp Count.csv")
inappres <- subset(inappres, select = c("Study", "Proportion.of.inapplicable", "Inapplicable"))
Resolution.morphy <- c(propnodemorphyAH2006,
                       propnodemorphyA2019,
                       propnodemorphyCPRS2019,
                       propnodemorphyC2013,
                       propnodemorphyC2020,
                       propnodemorphyD2016,
                       propnodemorphydL2020,
                       propnodemorphyE2014,
                       propnodemorphyFM2012,
                       propnodemorphyT2003)
inappres$Resolution.morphy <- Resolution.morphy

# when you have Anagallis data, colour by programme - may need to pivot_longer data frame

### proportion of inapp and resolution
pdf("Figures/is there a correlation between inapp proportion and resolution? morphy.pdf", 8, 5.83)
ggplot(inappres, aes(Proportion.of.inapplicable, Resolution.morphy)) +
  geom_point(size = 2, colour = palette[2]) +
  geom_smooth(se = F, method = "lm", fullrange = T, size = 0.75, colour = palette[2]) +
  theme_pubr() +
  labs(x = "Proportion of inapplicable data", y = "Resolution")
dev.off()

## correlation test
# resolution is normally distributed - for tests check "is there a relationship between resolution and distance from paup"
shapiro.test(inappdist$Proportion.of.inapplicable)
# proportion of inapp data is not normally distributed => use non-parametric test
corpropinappxres <- cor.test(inappres$Proportion.of.inapplicable, inappres$Resolution.morphy,
                             method = "spearman")
# rho=0.5151515, p=0.1328 => non-significant positive correlation


### inapp count and resolution
pdf("Figures/is there a correlation between inapp count and resolution? morphy.pdf", 8, 5.83)
ggplot(inappres, aes(Inapplicable, Resolution.morphy)) +
  geom_point(size = 2, colour = palette[2]) +
  geom_smooth(se = F, method = "lm", fullrange = T, size = 0.75, colour = palette[2]) +
  theme_pubr() +
  labs(x = "Inapplicable count", y = "Resolution")
dev.off()

## correlation test
shapiro.test(inappres$Inapplicable)
# inapplicable count is not normally distributed => non-parametric test
corinappcountxres <- cor.test(inappres$Inapplicable, inappres$Resolution.morphy,
                              method = "spearman")
# rho=0.5393939, p=0.1133 => non-significant positive correlation


#### how many times are the trees generated by morphy different to those generated by paup? ####
# tree contradiction
yestc <- 11
notc <- 1

# RF
yesrf <- 11
norf <- 1

# quartet
yesq <- 11
noq <- 1

yes <- c(yestc, yesrf, yesq)
no <- c(notc, norf, noq)
difftrees <- as.data.frame(cbind(yes, no))
rownames(difftrees) <- c("Tree contradiction", "Robinson-Foulds distance", "Quartet distance")
difftrees <- setDT(difftrees, keep.rownames = T)
difftrees <- pivot_longer(difftrees, cols = c("yes", "no"),
                          names_to = "Significant difference", values_to = "Count")

pdf("Figures/how many times is morphy significantly different to paup.pdf", 8, 5.83)
ggplot(difftrees) +
  geom_col(aes(rn, Count, fill = `Significant difference`), position = "dodge") +
  labs(x = "Tree difference metric", y = "Count") +
  theme_pubr() +
  scale_y_continuous(n.breaks = 6) +
  fill_palette(palette = c("firebrick", "steelblue"))
dev.off()  
  
