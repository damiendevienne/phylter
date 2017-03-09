library(DistatisR)
library(phangorn)
source("/media/aurore/KINGSTON/stageLBBE/pmcoa.R")
source("/home/aurore/Documents/Phylter/Fonctions1.R")
source("/home/aurore/Documents/Phylter/Simulation.R")

setwd(dir="/home/aurore/Documents/Phylter/trees/rose/")

#write.tree(rtree(20),"arbreHGT.tree") 

ListOut = SimOutliersHGT(nbsp = 50, nbgn = 10, outgn= 0, outsp = 3, sp=1)
RES <-Phylter(ListOut, distance="patristic", k=1.5, thres=0.5, quiet=TRUE)
plot.2WR(RES$Complete$mat2WR)
#plot(read.tree("arbre.tree"))
#RES <-pMCOA.complete(ListOut, distance="patristic", k=1.5, thres=0.5, quiet=TRUE)
RES$Complete$outgn
RES$Complete$outsp
RES$CellByCell$outcell


ListOut2 = SimOutliersLg(nbsp = 30, nbgn = 10, outgn= 1, outsp = 1, sp=1)
RES <-Phylter(ListOut2, distance="patristic", k=4, thres=0.5, quiet=TRUE)
plot.2WR(RES$Complete$mat2WR)
#RES <-pMCOA.complete(ListOut2, distance="patristic", k=1.5, thres=0.5, quiet=TRUE)
RES$Complete$outgn
RES$Complete$outsp
RES$CellByCell$outcell