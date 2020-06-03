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



### tree contradiction
tc <- sapply(treeslist, function(x)
  sapply(treeslist, treeContradiction, tree1=x))

colnames(tc) <- rownames(tc) <- c(paste("Paup", 1:17, sep = ""), 
                                  paste("Morphy", 1, sep = ""),
                                  paste("Anagallis", 1:165, sep = ""))
tc <- as.data.frame(tc)

# paup-paup
pptc <- tc[1:17, 1:17]
pptc <- as.matrix(pptc)
pptc <- as.vector(pptc)
pptcquartiles <- quantile(pptc)
pptcmedian <- pptcquartiles[3]

# paup-morphy
pmtc <- tc[18, 1:17]
pmtc <- as.matrix(pmtc)
pmtc <- as.numeric(pmtc)

# paup-anagallis
patc <- tc[1:17, 19:183]
patc <- as.matrix(patc)
patc <- as.vector(patc)

# morphy-anagallis
matc <- tc[18, 19:183]
matc <- as.matrix(matc)
matc <- as.vector(matc)


tclist <- list(pmtc, patc, matc)
tcquartiles <- as.data.frame(sapply(tclist, quantile, na.rm = T))
colnames(tcquartiles) <- c("PAUP* vs Morphy", "PAUP* vs Anagallis", "Morphy vs Anagallis")
tcmedian <- tcquartiles[3,]
tcLQ <- tcquartiles[2,]
tcUQ <- tcquartiles[4,]



### RF
rf <- sapply(treeslist, function(x)
  sapply(treeslist, RF.dist, tree1=x))

colnames(rf) <- rownames(rf) <- c(paste("Paup", 1:17, sep = ""), 
                                  paste("Morphy", 1, sep = ""),
                                  paste("Anagallis", 1:165, sep = ""))
rf <- as.data.frame(rf)

# paup-paup
pprf <- rf[1:17, 1:17]
pprf <- as.matrix(pprf)
pprf <- as.vector(pprf)
pprfquartiles <- quantile(pprf)
pprfmedian <- pprfquartiles[3]

# paup morphy
pmrf <- rf[18, 1:17]
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


### scale
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


### Quartets
qu <- sapply(treeslist, function(x)
  sapply(treeslist, QuartetPoints, cf=x))
qu1 <- as.data.frame(qu)
quartets <- c("Unresolved", "Contradicted", "Consistent")
quartets <- rep(quartets, length(qu1))
qu2 <- cbind(quartets, qu1)
qu3 <- subset(qu2, quartets == "Contradicted")
qu3 <- qu3[,2:ncol(qu3)]
colnames(qu3) <- rownames(qu3) <- c(paste("Paup", 1:17, sep = ""), 
                                  paste("Morphy", 1, sep = ""),
                                  paste("Anagallis", 1:165, sep = ""))
qu3 <- as.data.frame(qu3)
qu4 <- data.matrix(qu3)

# paup-paup
ppqu <- qu4[1:17, 1:17]
ppqu <- as.matrix(ppqu)
ppqu <- as.vector(ppqu)
ppququartiles <- quantile(ppqu)
ppqumedian <- ppququartiles[3]

# paup-morphy
pmqu <- qu4[18, 1:17]
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


### scale
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


#### MDS ####
# use raw distances
library(vegan)
mdsrf1 <- metaMDS(rf, distance = "euclidean", trymax = 999)
stressplot(mdsrf1)
plot(mdsrf1$points[,1],
     mdsrf1$points[,2],
     col = c("firebrick", "steelblue", "purple")[rf1$Programme])

# programme coloumn
paupnames <- c("PAUP*")
paupnames <- rep(paupnames, length(pauptrees))
morphynames <- c("Morphy")
morphynames <- rep(morphynames, length(morphytrees))
ananames <- c("Anagallis")
ananames <- rep(ananames, length(anatrees))
Programme <- c(paupnames, morphynames, ananames)

mdsrfpoints <- as.data.frame(mdsrf1$points)
mdsrfpoints <- cbind(mdsrfpoints, Programme)

palette <- wesanderson::wes_palette(name = "Chevalier1")
jitter <- position_jitter(width = 0.6, height = 0.6)

# Scatter plot with marginal density distribution plot
a <- ggplot(mdsrfpoints) +
  geom_point(aes(MDS1, MDS2, colour = Programme), position = jitter, size = 1.5, alpha = 1) +
  theme_pubr() +
  color_palette(palette = palette)


b <- ggplot(mdsrfpoints) +
  geom_density(aes(MDS1, fill = Programme), alpha = 0.5) +
  theme_pubr() +
  fill_palette(palette = palette)

c <- ggplot(mdsrfpoints) +
  geom_density(aes(MDS2, fill = Programme), alpha = 0.5) +
  theme_pubr() +
  fill_palette(palette = palette)

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
        axis.ticks = element_blank()
)

gridExtra::grid.arrange(b, blankPlot, a, c,
                        ncol = 2, nrow = 2, widths = c(4,1.4), heights = c(1.4,4))


# stat_ellipse(aes(MDS1, MDS2, colour = Programme))
ggplot(mdsrfpoints) +
  geom_point(aes(MDS1, MDS2, colour = Programme)) +
  stat_ellipse(aes(MDS1, MDS2, colour = Programme))


