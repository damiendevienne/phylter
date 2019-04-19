require(ape)
require(RSpectra)

rm(list=ls())
source("R/trees2matrices.R")
source("R/rename.genes.R")
source("R/impMean.R")
source("R/DistatisFast.R")
source("R/Dist2WR.R")
source("R/normalize.R")
source("R/detect.outliers.R")
source("R/phylter.R")
source("R/simulate.trees.R")

# #one dataset
# trees<-read.tree("/home/ddevienne/Documents/WORK/Phylter2018/supertree/mamma.tre")
# class(trees)<-"multiPhylo"
# gene.names<-NULL
# trees <- rename.genes(trees, gene.names = gene.names)



# trees<-read.tree("ArchaeaTrees.tre")
# class(trees)<-"multiPhylo"
# gene.names<-NULL
# trees <- rename.genes(trees, gene.names = gene.names)
