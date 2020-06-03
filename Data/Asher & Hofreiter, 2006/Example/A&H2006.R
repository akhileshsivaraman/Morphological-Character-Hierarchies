##### Asher & Hofreiter, 2006

setwd("~/Documents/Imperial/Morphological Character Hierarchies/Data/Asher & Hofreiter, 2006")
library(ggtree)
library(distory)
library(treeio)
library(ape)
library(phangorn)
library(tibble)
library(ggplot2)
library(phytools)
library(paleotree)
library(knitr)


paup <- read.nexus("paup A&H2006.tre")
paup <- consensus(paup, p=0.5)
is.rooted(paup)
paup <- root(paup, "Didelphis")

morphy <- read.nexus("Morphy A&H.tre")
is.rooted(morphy)

ana <- read.tree(text = "(Didelphis,((Canis, Procavia),(((((Crocidura, Blarina),Sorex),((Chrysochloris,((((Microgale, Oryzorictes),((Potamogale, Micropotamogale),Limnogale)),(Erythrozootes,((Geogale, Parageogale),Protenrec))),((Hemicentetes, Tenrec),(Echinops, Setifer)))),Erinaceus)),Elephantulus),Orycteropus)));")
is.rooted(ana)

par(mar = c(0,0,0,0))
plot.phylo(paup, type = "cladogram")
plot.phylo(morphy, type = "cladogram")
plot.phylo(ana, type = "cladogram")


morphy <- rotateNodes(morphy, "all")
ana <- rotateNodes(ana, "all")
plot.phylo(paup, type = "cladogram")
plot.phylo(morphy, type = "cladogram")
plot.phylo(ana, type = "cladogram")

##### graphical comparison
comparePhylo(paup, morphy)
comparePhylo(paup, ana)
comparePhylo(morphy, ana)
# can I put them into a little table?

par(mar = c(0,0,0,0))
phylo.diff(paup, morphy, type = "cladogram", edge.width = 2.5, label.offset = 0.5)
phylo.diff(paup, ana, type = "cladogram", edge.width = 2.5, label.offset = 0.5)
phylo.diff(morphy, ana, type = "cladogram", edge.width = 2.5, label.offset = 0.5)


par(mfrow = c(1,1))
compare.chronograms(paup, morphy, type = "cladogram")
par(mfrow = c(1,2))
plot.phylo(paup)
plot.phylo(morphy)

##### testing contradictions
treelist <- c(morphy,paup,ana)

dissimilarity <- sapply(treelist, function(x)
  sapply(treelist, treeContradiction, tree1=x))

randomtreemean <- function(ntaxa, replicates = 100){
  return(mean(unlist(mapply(treeContradiction, rmtree(replicates, ntaxa), 
                            rmtree(replicates, ntaxa), SIMPLIFY = F))))
}
set.seed(111)

tmeanA <- randomtreemean(Ntip(treelist[[1]]))

correction <- function(dissimilarity, tmean){
  return((tmean - dissimilarity)/tmean)
}

similarity <- apply(dissimilarity, c(1,2), correction, tmeanA)

colnames(similarity) <- rownames(similarity) <- c("Morphy", "Paup", "Anagallis")
table <- kable(similarity, caption = "similarity")
kableExtra::save_kable(table, file = "A&H2006.txt")


##### Robsinson-Foulds
multiRF(treelist)
is.binary(paup)
is.binary(morphy)
is.binary(ana)

# randomly resolve polytomies in paup
paup1 <- multi2di(paup)
trees <- c(paup1, morphy, ana)
multiRF(trees)

treedist(paup1, morphy)
treedist(morphy, ana)
