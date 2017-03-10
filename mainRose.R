library(DistatisR)
library(phangorn)
source("/media/aurore/KINGSTON/stageLBBE/pmcoa.R")
source("/home/aurore/Documents/Phylter/Fonctions1.R")
source("/home/aurore/Documents/Phylter/Simulation.R")

setwd(dir="/home/aurore/Documents/Phylter/trees/rose/")

#write.tree(rtree(20),"arbreHGT.tree") 

ListOut = SimOutliersHGT(nbsp = 50, nbgn = 10, outgn= 1, outsp = 2, sp=1)
RESTree <-Phylter(ListOut$ListTrees, distance="nodal", k=2, thres=0.5, quiet=TRUE)
RESSim <-Phylter(ListOut$ListSim, distance="nodal", k=2, thres=0.5, quiet=TRUE)

plot.2WR(RESTree$Complete$mat2WR)
plot.2WR(RESSim$Complete$mat2WR)

RESTree$Complete$outgn
RESTree$Complete$outsp

RESSim$Complete$outgn
RESSim$Complete$outsp


#plot(read.tree("arbre.tree"))
#tree=read.tree("arbre.tree")
#RES <-pMCOA.complete(ListOut, distance="nodal", k=1.5, thres=0.5, quiet=TRUE)
#plot.2WR(RES$step1$mat2WR)

ListOut2 = SimOutliersLg(nbsp = 10, nbgn = 10, outgn= 1, outsp = 1, sp=1)

RESTree <-Phylter(ListOut2$ListTree, distance="patristic", k=2, thres=0.5, quiet=TRUE)
RESSim <-Phylter(ListOut2$ListSim, distance="patristic", k=2, thres=0.5, quiet=TRUE)

plot.2WR(RESTree$Complete$mat2WR)
plot.2WR(RESSim$Complete$mat2WR)

RES$Complete$outgn
RES$Complete$outsp
RES$CellByCell$outcell