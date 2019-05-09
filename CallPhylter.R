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
source("R/summary.phylter.R")
source("R/plot2WR.R")


genes<-readLines("ALLNAMES0.5")
trees<-list()
for (i in 1:length(genes)) {
	trees[[i]]<-read.tree(genes[i])
}

# #one dataset
trees<-read.tree("/home/ddevienne/Documents/WORK/Phylter2018/supertree/mamma.tre")
class(trees)<-"multiPhylo"
gene.names<-NULL
trees <- rename.genes(trees, gene.names = gene.names)



trees<-read.tree("ArchaeaTrees.tre")
class(trees)<-"multiPhylo"
gene.names<-NULL
trees <- rename.genes(trees, gene.names = gene.names)


