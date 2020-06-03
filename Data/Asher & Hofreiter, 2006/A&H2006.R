##### Asher & Hofreiter, 2006 #####

setwd("~/Documents/Imperial/Morphological Character Hierarchies/Data/Asher & Hofreiter, 2006")
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
pauptrees <- read.nexus("Paup-A&H2006.tre")
paup <- consensus(pauptrees)
paup <- unroot(paup)
pauptrees <- as.multiPhylo(pauptrees)
pauptrees <- unroot.multiPhylo(pauptrees)

# morphy
morphytrees <- read.nexus("Morphy-A&H2006.tre")
morphy <- consensus(morphytrees)
morphy <- unroot(morphy)
morphytrees <- as.multiPhylo(morphytrees)
morphytrees <- unroot.multiPhylo(morphytrees)

# anagallis
anatrees <- read.nexus("Anatrees-A&H2006.txt.nex")
ana <- consensus(anatrees)
ana <- unroot(ana)
anatrees <- list.filter(anatrees, Nnode == 22)
anatrees <- as.multiPhylo(anatrees)
anatrees <- unroot.multiPhylo(anatrees)


contreelist <- c(paup, morphy, ana)
treeslist <- c(pauptrees, morphytrees, anatrees)


#### plot contrees ####
par(mar=c(0,0,0,0))
par(mfrow=c(1,1))
plot.phylo(paup, type = "cladogram")
plot.phylo(morphy, type = "cladogram")
plot.phylo(ana, type = "cladogram")



#### basic comparison of contrees ####
comparePhylo(paup, morphy)
comparePhylo(paup, ana)
comparePhylo(morphy, ana)

# visual comparison
par(mar=c(0,0,0,0))
pdf("Figures/phylodiff paup vs morphy.pdf", 8.27, 5.83)
phylo.diff(paup, morphy, type = "cladogram", edge.width = 1.5, label.offset = 0.5, no.margin = T)
dev.off()
pdf("Figures/phylodiff paup vs ana.pdf", 8.27, 5.83)
phylo.diff(paup, ana, type = "cladogram", edge.width = 1.5, label.offset = 0.5, no.margin = T)
dev.off()
pdf("Figures/phylodiff morphy vs ana.pdf", 8.27, 5.83)
phylo.diff(morphy, ana, type = "cladogram", edge.width = 1.5, label.offset = 0.5, no.margin = T)
dev.off()


##### Metrics #####
#### RF ####
rf <- sapply(treeslist, function(x)
  sapply(treeslist, RF.dist, tree1=x))

colnames(rf) <- rownames(rf) <- c(paste("Paup", 1:length(pauptrees), sep = ""), 
                                  paste("Morphy", 1:length(morphytrees), sep = ""),
                                  paste("Anagallis", 1:length(anatrees), sep = ""))
rf <- as.data.frame(rf)

# paup-paup
pprf <- rf[1:17, 1:17]
pprf <- as.matrix(pprf)
pprf <- as.vector(pprf)
pprfquartiles <- quantile(pprf)
pprfmedian <- pprfquartiles[3]

# morphy-morphy
mmrf <- rf[18, 18]
mmrf <- as.matrix(mmrf)
mmrf <- as.vector(mmrf)
mmrfquartiles <- quantile(mmrf)
mmrfmedian <- mmrfquartiles[3]

# ana-ana
aarf <- rf[19:183, 19:183]
aarf <- as.matrix(aarf)
aarf <- as.vector(aarf)
aarfquartiles <- quantile(aarf)
aarfmedian <- aarfquartiles[3]

# paup morphy
pmrf <- rf[1:17, 18]
pmrf <- as.matrix(pmrf)
pmrf <- as.vector(pmrf)

# paup anagallis
parf <- rf[1:17, 19:183]
parf <- as.matrix(parf)
parf <- as.vector(parf)

# morphy anagallis
marf <- rf[18, 19:183]
marf <- as.matrix(marf)
marf <- as.vector(marf)

# median
rflist <- list(pmrf, parf, marf)
rfquartiles <- as.data.frame(sapply(rflist, quantile, na.rm = T))
colnames(rfquartiles) <- c("PAUP* vs Morphy", "PAUP* vs Anagallis", "Morphy vs Anagallis")
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
                                    paste("Morphy", 1:length(morphytrees), sep = ""),
                                    paste("Anagallis", 1:length(anatrees), sep = ""))
qu3 <- as.data.frame(qu3)
qu4 <- data.matrix(qu3)

# paup-paup
ppqu <- qu4[1:17, 1:17]
ppqu <- as.matrix(ppqu)
ppqu <- as.vector(ppqu)
ppququartiles <- quantile(ppqu)
ppqumedian <- ppququartiles[3]

# morphy-morphy
mmqu <- qu4[18, 18]
mmqu <- as.matrix(mmqu)
mmqu <- as.vector(mmqu)
mmququartiles <- quantile(mmqu)
mmqumedian <- mmququartiles[3]

# ana-ana
aaqu <- qu4[19:183, 19:183]
aaqu <- as.matrix(aaqu)
aaqu <- as.vector(aaqu)
aaququartiles <- quantile(aaqu)
aaqumedian <- aaququartiles[3]

# paup-morphy
pmqu <- qu4[1:17, 18]
pmqu <- as.matrix(pmqu)
pmqu <- as.vector(pmqu)
pmququartiles <- quantile(pmqu)


# paup-anagallis
paqu <- qu4[1:17, 19:183]
paqu <- as.matrix(paqu)
paqu <- as.vector(paqu)
paququartiles <- quantile(paqu)

# morphy-anagallis
maqu <- qu4[18, 19:183]
maqu <- as.matrix(maqu)
maqu <- as.vector(maqu)
maququartiles <- quantile(maqu)

# median and quartiles
qulist <- list(pmqu, paqu, maqu)
ququartiles <- as.data.frame(sapply(qulist, quantile, na.rm = T))
colnames(ququartiles) <- c("PAUP* vs Morphy", "PAUP* vs Anagallis", "Morphy vs Anagallis")
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
rownames(rawdata) <- c("Paup vs Morphy", "Paup vs Anagallis", "Morphy vs Anagallis")
rawdata$`Robinson-Foulds Distance UQ` <- t(rfUQ)
rawdata$`Robinson-Foulds Distance LQ` <- t(rfLQ)
rawdata$`Quartet Distances UQ` <- t(quUQ)
rawdata$`Quartet Distances LQ` <- t(quLQ)
write.csv(rawdata, file = "rawdataAH2006.csv")

# normalised data
normdata <- as.data.frame(t(rbind(normalisedRFmedian, normalisedQmedian)))
colnames(normdata) <- c("Robinson-Foulds Distance", "Quartet Distances")
rownames(normdata) <- c("Paup vs Morphy", "Paup vs Anagallis", "Morphy vs Anagallis")
normdata$`Robinson-Foulds Distance UQ` <- t(normalisedRFUQ)
normdata$`Robinson-Foulds Distance LQ` <- t(normalisedRFLQ)
normdata$`Quartet Distances UQ` <- t(normalisedQUQ)
normdata$`Quartet Distances LQ` <- t(normalisedQLQ)
write.csv(normdata, file = "normdataAH2006.csv")



##### Stats tests ##### redo 
#### Bhattacharyya ####
# measuring similarity among distributions
# probability of overlap between distributions
# is the distribution of this pairwise comparison similar to the distribution of that pairwise comparison?
# <0.05 significantly different, >0.95 significantly similar
# probability that these trees are as different from random trees as those trees


# RF - normalise distributions to use Bhattacharyya
normpprf <- as.numeric(pprf/randoRFmedian)
normpmrf <- as.numeric(pmrf/randoRFmedian)
normparf <- as.numeric(parf/randoRFmedian)
normmarf <- as.numeric(marf/randoRFmedian)

bhatt.coeff(normpprf, normpmrf) # 0.4459379
bhatt.coeff(normpprf, normparf) # 0.9327063
bhatt.coeff(normpmrf, normparf) # 0.5952563
bhatt.coeff(normpmrf, normmarf) # 0.742394
bhatt.coeff(normparf, normmarf) # 0.6750986


# Quartets
normppqu <- as.numeric(ppqu/randoQmedian)
normpmqu <- as.numeric(pmqu/randoQmedian)
normpaqu <- as.numeric(paqu/randoQmedian)
normmaqu <- as.numeric(maqu/randoQmedian)

bhatt.coeff(normppqu, normpmqu) # 0.347582
bhatt.coeff(normppqu, normpaqu) # 0.943885
bhatt.coeff(normpmqu, normpaqu) # 0.5789239
bhatt.coeff(normpmqu, normmaqu) # 0.8234817
bhatt.coeff(normpaqu, normmaqu) # 0.549577


#### t-test/Wilcoxon ####
## which programme has sampled the largest area of the tree space?

### RF
# paup-morphy
wilcox.test(pprf, mmrf) # W = 280.5, p-value = 0.1012

# paup-ana
wilcox.test(pprf, aarf) # W = 3883468, p-value = 0.7051

# morphy-ana
wilcox.test(mmrf, aarf) # W = 724.5, p-value = 0.09914


### Quartet
# paup-morphy
wilcox.test(ppqu, mmqu) # W = 280.5, p-value = 0.1012

# paup-ana
wilcox.test(ppqu, aaqu) # W = 3473860, p-value = 0.0006118
# aaqu > ppqu

# morphy-ana
wilcox.test(mmqu, aaqu) # W = 724.5, p-value = 0.101


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
ananames <- c("Anagallis")
ananames <- rep(ananames, length(anatrees))
Programme <- c(paupnames, morphynames, ananames)


### RF
mdsrf <- metaMDS(rf, distance = "euclidean", trymax = 999)
stressplot(mdsrf)
mdsrfpoints <- as.data.frame(mdsrf$points)
mdsrfpoints <- cbind(mdsrfpoints, Programme)


# Scatter plot with marginal density distribution plot
jitter <- position_jitter(width = 1.5, height = 1)
mdsrfplot1 <- ggplot(mdsrfpoints) +
  geom_point(aes(MDS1, MDS2, colour = Programme), position = jitter, size = 2.5, alpha = 1) +
  theme_pubr() + 
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.title = element_blank(), 
        axis.text = element_blank(), legend.position = "none", panel.border = element_rect(colour = "grey", size = 0.5, fill = NA)) +
  color_palette(palette = c(palette[3], palette[2], palette[1]))

mdsrfplot2 <- ggplot(mdsrfpoints) +
  geom_density(aes(MDS1, fill = Programme), alpha = 0.4) +
  theme_pubr() +
  theme(axis.title = element_text(size = 11)) +
  fill_palette(palette = c(palette[3], palette[2], palette[1])) +
  scale_y_continuous(breaks = c(0,0.06))

mdsrfplot3 <- ggplot(mdsrfpoints) +
  geom_density(aes(y = MDS2, fill = Programme), alpha = 0.4) +
  theme_pubr() +
  theme(legend.position = "none", axis.title = element_text(size = 11)) +
  fill_palette(palette = c(palette[3], palette[2], palette[1]))

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

pdf("Figures/RF MDS with density distribution AH2006.pdf", 8, 5.3)
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


# Scatter plot with marginal density distribution plot
jitter <- position_jitter(width = 200, height = 200)
mdsquplot1 <- ggplot(mdsqupoints) +
  geom_point(aes(MDS1, MDS2, colour = Programme), position = jitter, size = 2.5, alpha = 1) +
  theme_pubr() + 
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.title = element_blank(), 
        axis.text = element_blank(), legend.position = "none", panel.border = element_rect(colour = "grey", size = 0.5, fill = NA)) +
  color_palette(palette = c(palette[3], palette[2], palette[1]))

mdsquplot2 <- ggplot(mdsqupoints) +
  geom_density(aes(MDS1, fill = Programme), alpha = 0.4) +
  theme_pubr() +
  theme(axis.title = element_text(size = 11)) +
  fill_palette(palette = c(palette[3], palette[2], palette[1])) +
  scale_y_continuous(breaks = c(0, 0.0005))

mdsquplot3 <- ggplot(mdsqupoints) +
  geom_density(aes(y = MDS2, fill = Programme), alpha = 0.4) +
  theme_pubr() +
  theme(legend.position = "none", axis.title = element_text(size = 11)) +
  fill_palette(palette = c(palette[3], palette[2], palette[1]))


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

pdf("Figures/Quartet MDS with density distribution AH2006.pdf", 8, 5.3)
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


