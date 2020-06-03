##### Lihoreau et al, 2015 #####

setwd("~/Documents/Imperial/Morphological Character Hierarchies/Data/Lihoreau et al, 2015/")
library(ape)
library(distory)
library(phangorn)
library(phytools)
library(paleotree)
library(dplyr)
library(Quartet)
library(dispRity)
library(rlist)
library(ggpubr)


# paup
pauptrees <- read.nexus("Pauptrees-L2015.tre")
paup <- read.nexus("Paupcontree-L2015.tre")
paup <- unroot(paup)
pauptrees <- as.multiPhylo(pauptrees)
pauptrees <- unroot.multiPhylo(pauptrees)

# morphy
morphytrees <- read.nexus("Morphytrees-L2015.tre")
morphy <- read.nexus("Morphycontree-L2015.tre")
morphy <- unroot(morphy)
morphytrees <- as.multiPhylo(morphytrees)
morphytrees <- unroot.multiPhylo(morphytrees)


contreelist <- c(paup, morphy)
treeslist <- c(pauptrees, morphytrees)


#### plot contrees ####
par(mar=c(0,0,0,0))
par(mfrow=c(1,1))
plot.phylo(paup, type = "cladogram")
plot.phylo(morphy, type = "cladogram")



#### basic comparison of contrees ####
comparePhylo(paup, morphy)

# visual comparison
par(mar=c(0,0,0,0))
pdf("Figures/phylodiff paup vs morphy.pdf", 8.27, 5.83)
phylo.diff(paup, morphy, type = "cladogram", edge.width = 1.5, label.offset = 0.5, no.margin = T)
dev.off()


##### Metrics #####
#### RF ####
rf <- sapply(treeslist, function(x)
  sapply(treeslist, RF.dist, tree1=x))

colnames(rf) <- rownames(rf) <- c(paste("Paup", 1:length(pauptrees), sep = ""), 
                                  paste("Morphy", 1:length(morphytrees), sep = ""))
rf <- as.data.frame(rf)

# paup-paup
pprf <- rf[1:2, 1:2]
pprf <- as.matrix(pprf)
pprf <- as.vector(pprf)
pprfquartiles <- quantile(pprf)
pprfmedian <- pprfquartiles[3]

# morphy-morphy
mmrf <- rf[3:6, 3:6]
mmrf <- as.matrix(mmrf)
mmrf <- as.vector(mmrf)
mmrfquartiles <- quantile(mmrf)
mmrfmedian <- mmrfquartiles[3]


# paup morphy
pmrf <- rf[1:2, 3:6]
pmrf <- as.matrix(pmrf)
pmrf <- as.vector(pmrf)


# median
rflist <- list(pmrf)
rfquartiles <- as.data.frame(sapply(rflist, quantile, na.rm = T))
colnames(rfquartiles) <- c("PAUP* vs Morphy")
rfmedian <- rfquartiles[3,]
rfLQ <- rfquartiles[2,]
rfUQ <- rfquartiles[4,]


### normalise
random.comparisons.quantile <- function(number.tips, number.trees, difference.function) {
  trees1 <- rmtree(number.trees, number.tips, br = NULL, rooted = F)
  trees2 <- rmtree(number.trees, number.tips, br = NULL, rooted = F)
  result <- quantile(unlist(mapply(difference.function, trees1, trees2, SIMPLIFY = F)))
  return(result)
}
set.seed(12345)
randoRFquantile <- random.comparisons.quantile(Ntip(paup), 100, RF.dist)
randoRFmedian <- randoRFquantile[3]
randoRFLQ <- randoRFquantile[2]
randoRFUQ <- randoRFquantile[4]


normalisedRFmedian <- rfmedian/randoRFmedian
normalisedRFLQ <- rfLQ/randoRFLQ
normalisedRFUQ <- rfUQ/randoRFUQ



#### Quartet ####
qu <- sapply(treeslist, function(x)
  sapply(treeslist, QuartetPoints, cf=x))
qu1 <- as.data.frame(qu)
quartets <- c("Unresolved", "Contradicted", "Consistent")
quartets <- rep(quartets, length(qu1))
qu2 <- cbind(quartets, qu1)
qu3 <- subset(qu2, quartets == "Contradicted")
qu3 <- qu3[,2:ncol(qu3)]
colnames(qu3) <- rownames(qu3) <- c(paste("Paup", 1:length(pauptrees), sep = ""), 
                                    paste("Morphy", 1:length(morphytrees), sep = ""))
qu3 <- as.data.frame(qu3)
qu4 <- data.matrix(qu3)

# paup-paup
ppqu <- qu4[1:2, 1:2]
ppqu <- as.matrix(ppqu)
ppqu <- as.vector(ppqu)
ppququartiles <- quantile(ppqu)
ppqumedian <- ppququartiles[3]

# morphy-morphy
mmqu <- qu4[3:6, 3:6]
mmqu <- as.matrix(mmqu)
mmqu <- as.vector(mmqu)
mmququartiles <- quantile(mmqu)
mmqumedian <- mmququartiles[3]

# paup-morphy
pmqu <- qu4[1:2, 3:6]
pmqu <- as.matrix(pmqu)
pmqu <- as.vector(pmqu)
pmququartiles <- quantile(pmqu)


# median and quartiles
qulist <- list(pmqu)
ququartiles <- as.data.frame(sapply(qulist, quantile, na.rm = T))
colnames(ququartiles) <- c("PAUP* vs Morphy")
qumedian <- ququartiles[3,]
quLQ <- ququartiles[2,]
quUQ <- ququartiles[4,]


### normalise
random.Q.quantile <- function(number.tips, number.trees, difference.function) {
  trees1 <- rmtree(number.trees, number.tips, br = NULL, rooted = F)
  trees2 <- rmtree(number.trees, number.tips, br = NULL, rooted = F)
  result <- mapply(difference.function, trees1, trees2, SIMPLIFY = F)
  result <- bind_rows(result)
  result <- sapply(result, quantile)
  return(result)
}
set.seed(12345)
randoQquantile <- random.Q.quantile(Ntip(paup), 100, QuartetPoints)
randoQquantile <- randoQquantile[,2]
randoQmedian <- randoQquantile[3]
randoQLQ <- randoQquantile[2]
randoQUQ <- randoQquantile[4]

normalisedQmedian <- qumedian/randoQmedian
normalisedQLQ <- quLQ/randoQLQ
normalisedQUQ <- quUQ/randoQUQ



#### save data ####
# raw data - note: Tree Contradiction is already normalised
rawdata <- as.data.frame(t(rbind(rfmedian, qumedian)))
colnames(rawdata) <- c("Robinson-Foulds Distance", "Quartet Distances")
rownames(rawdata) <- c("Paup vs Morphy")
rawdata$`Robinson-Foulds Distance UQ` <- t(rfUQ)
rawdata$`Robinson-Foulds Distance LQ` <- t(rfLQ)
rawdata$`Quartet Distances UQ` <- t(quUQ)
rawdata$`Quartet Distances LQ` <- t(quLQ)
write.csv(rawdata, file = "rawdataL2015.csv")

# normalised data
normdata <- as.data.frame(t(rbind(normalisedRFmedian, normalisedQmedian)))
colnames(normdata) <- c("Robinson-Foulds Distance", "Quartet Distances")
rownames(normdata) <- c("Paup vs Morphy")
normdata$`Robinson-Foulds Distance UQ` <- t(normalisedRFUQ)
normdata$`Robinson-Foulds Distance LQ` <- t(normalisedRFLQ)
normdata$`Quartet Distances UQ` <- t(normalisedQUQ)
normdata$`Quartet Distances LQ` <- t(normalisedQLQ)
write.csv(normdata, file = "normdataL2015.csv")



#### Bhattacharyya ####
# measuring similarity among distributions
# probability of overlap between distributions
# is the distribution of this pairwise comparison similar to the distribution of that pairwise comparison?
# <0.05 significantly different, >0.95 significantly similar
# probability that these trees are as different from random trees as those trees


# RF - normalise distributions to use Bhattacharyya
normpprf <- as.numeric(pprf/randoRFmedian)
normpmrf <- as.numeric(pmrf/randoRFmedian)

bhatt.coeff(normpprf, normpmrf) # 0


# Quartets
normppqu <- as.numeric(ppqu/randoQmedian)
normpmqu <- as.numeric(pmqu/randoQmedian)

bhatt.coeff(normppqu, normpmqu) # 0


#### t-test/Wilcoxon ####
## which programme has sampled the largest area of the tree space?

### RF
# paup-morphy
wilcox.test(pprf, mmrf) # W = 16, p-value = 0.1293
# neither larger

### Quartet
# paup-morphy
wilcox.test(ppqu, mmqu) # W = 12, p-value = 0.05875
# neither larger


##### Visualisation #####
palette <- wesanderson::wes_palette(name = "Chevalier1")
palette <- c(palette, "#25291C", "#E6E49F")
# PAUP* = palette[1]
# Morphy = palette[2]
# Anagallis = palette[3]
# PAUP* vs Morphy = palette[4]
# PAUP* vs Anagallis = palette[5]
# Morphy vs Anagallis = palette[6]
comparepalette <- c(palette[6], palette[5], palette[4])

#### Median and IQR ####
library(data.table)
library(gridExtra)
library(cowplot)
par(mar=c(0,0,0,0))

rawdata <- setDT(rawdata, keep.rownames = T)
# RF plot
rfplot <- ggplot(rawdata, aes(rn, rawdata$`Robinson-Foulds Distance`)) + theme_pubr() +
  geom_pointrange(aes(x = rn, ymin = rawdata$`Robinson-Foulds Distance LQ`, ymax = rawdata$`Robinson-Foulds Distance UQ`,
                      colour = rn)) +
  labs(x = "Pairwise comparison",
       y = "Robinson-Foulds 
distance") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title = element_text(size = 12)) +
  color_palette(palette = comparepalette)
# Quartet
quartplot <- ggplot(rawdata, aes(rn, rawdata$`Quartet Distances`)) + theme_pubr() +
  geom_pointrange(aes(x = rn, ymin = rawdata$`Quartet Distances LQ`, ymax = rawdata$`Quartet Distances UQ`,
                      colour = rn)) +
  labs(x = "Pairwise comparison",
       y = "Quartet distance") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title = element_text(size = 12),
        legend.box = "horizontal", legend.position = "bottom", legend.title = element_blank()) +
  color_palette(palette = comparepalette)

averageplots <- plot_grid(rfplot + theme(legend.position = "none"),
                          quartplot + theme(legend.position = "none"),
                          labels = "AUTO", nrow = 1, ncol = 2)
legend <- get_legend(quartplot + 
                       theme(legend.position = "bottom") + 
                       guides(colour = guide_legend(nrow = 1)))

pdf("Figures/median.pdf", 8, 5.83)
plot_grid(averageplots, legend, ncol = 1, rel_heights = c(1,0.1))
dev.off()


#### MDS ####
library(vegan)
library(cowplot)
# programme column
paupnames <- c("PAUP*")
paupnames <- rep(paupnames, length(pauptrees))
morphynames <- c("Morphy")
morphynames <- rep(morphynames, length(morphytrees))
Programme <- c(paupnames, morphynames)


### RF
mdsrf <- metaMDS(rf, distance = "euclidean", trymax = 999)
stressplot(mdsrf)
mdsrfpoints <- as.data.frame(mdsrf$points)
mdsrfpoints <- cbind(mdsrfpoints, Programme)


# Scatter plot with marginal density distribution plot
mdsrfplot1 <- ggplot(mdsrfpoints) +
  geom_point(aes(MDS1, MDS2, colour = Programme), position = "jitter", size = 2.5, alpha = 1) +
  theme_pubr() + 
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.title = element_blank(), 
        axis.text = element_blank(), legend.position = "none", panel.border = element_rect(colour = "grey", size = 0.5, fill = NA)) +
  color_palette(palette = c(palette[2], palette[1]))

mdsrfplot2 <- ggplot(mdsrfpoints) +
  geom_density(aes(MDS1, fill = Programme), alpha = 0.4) +
  theme_pubr() +
  theme(axis.title = element_text(size = 11)) +
  fill_palette(palette = c(palette[2], palette[1])) +
  scale_y_continuous(breaks = c(0, 0.3))

mdsrfplot3 <- ggplot(mdsrfpoints) +
  geom_density(aes(y=MDS2, fill = Programme), alpha = 0.4) +
  theme_pubr() +
  theme(legend.position = "none", axis.title = element_text(size = 11)) +
  fill_palette(palette = c(palette[2], palette[1])) +
  scale_x_continuous(breaks = c(0, 9))

blankPlot <- ggplot()+geom_blank(aes(1,1))+
  theme(plot.background = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(), 
        axis.text.y = element_blank(),
        axis.ticks = element_blank())

pdf("Figures/RF MDS with density distribution L2015.pdf", 8, 5.3)
ggdraw() +
  draw_plot(mdsrfplot2, x = 0, y = 0.7, width = 0.7, height = 0.3) +
  draw_plot(mdsrfplot1, x = 0, y = 0, width = 0.7, height = 0.7) +
  draw_plot(mdsrfplot3, x = 0.7, y = 0, width = 0.3, height = 0.7) +
  draw_plot_label(label = "A", x = 0, y = 0.98) +
  draw_plot_label(label = "B", x = 0, y = 0.74) +
  draw_plot_label(label = "C", x = 0.7, y = 0.74)
dev.off()


### Quartet
mdsqu <- metaMDS(qu4, distance = "euclidean", trymax = 999)
stressplot(mdsqu)
mdsqupoints <- as.data.frame(mdsqu$points) 
mdsqupoints <- cbind(mdsqupoints, Programme)

jitter <- position_jitter(width = 0.6, height = 0.6)

# Scatter plot with marginal density distribution plot
mdsquplot1 <- ggplot(mdsqupoints) +
  geom_point(aes(MDS1, MDS2, colour = Programme), position = "jitter", size = 2.5, alpha = 1) +
  theme_pubr() + 
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.title = element_blank(), 
        axis.text = element_blank(), legend.position = "none", panel.border = element_rect(colour = "grey", size = 0.5, fill = NA)) +
  color_palette(palette = c(palette[2], palette[1]))

mdsquplot2 <- ggplot(mdsqupoints) +
  geom_density(aes(MDS1, fill = Programme), alpha = 0.4) +
  theme_pubr() +
  theme(axis.title = element_text(size = 11)) +
  fill_palette(palette = c(palette[2], palette[1])) +
  scale_y_continuous(breaks = c(0, 0.0075))

mdsquplot3 <- ggplot(mdsqupoints) +
  geom_density(aes(y=MDS2, fill = Programme), alpha = 0.4) +
  theme_pubr() +
  theme(legend.position = "none", axis.title = element_text(size = 11)) +
  fill_palette(palette = c(palette[2], palette[1])) +
  scale_x_continuous(breaks = c(0, 20))


blankPlot <- ggplot()+geom_blank(aes(1,1))+
  theme(plot.background = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(), 
        axis.text.y = element_blank(),
        axis.ticks = element_blank())

pdf("Figures/Quartet MDS with density distribution L2015.pdf", 8, 5.3)
ggdraw() +
  draw_plot(mdsquplot2, x = 0, y = 0.7, width = 0.7, height = 0.3) +
  draw_plot(mdsquplot1, x = 0, y = 0, width = 0.7, height = 0.7) +
  draw_plot(mdsquplot3, x = 0.7, y = 0, width = 0.3, height = 0.7) +
  draw_plot_label(label = "A", x = 0, y = 0.98) +
  draw_plot_label(label = "B", x = 0, y = 0.74) +
  draw_plot_label(label = "C", x = 0.7, y = 0.74)
dev.off()


#### Pair plot ####
par(mfrow = c(1,1))
par(mar=c(2.5,2.5,2.5,0.5))
pair.plot(normdata, what = 1, col = c("#90AFC5","#0076C9"), diag = 1, main = "Tree contradiction")
legend("topleft",
       legend = c("1: PAUP* vs Morphy", 
                  "2: PAUP* vs Anagallis", 
                  "3: Morphy vs Anagallis"),
       cex = 0.5,
       box.lty = 0)
pair.plot(normdata, what = 2, col = c("#90AFC5","#0076C9"), diag = 1, legend = T)
pair.plot(normdata, what = 3, col = c("#90AFC5","#0076C9"), diag = 1, legend = T)
# add labels - ask Thomas


