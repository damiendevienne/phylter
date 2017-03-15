library(DistatisR)
library(phangorn)
source("/media/aurore/KINGSTON/stageLBBE/pmcoa.R")
source("/home/aurore/Documents/Phylter/Fonctions1.R")
source("/home/aurore/Documents/Phylter/Simulation.R")

setwd(dir="/home/aurore/Documents/Phylter/trees/rose/")

ListOut = SimOutliersHGT(nbsp = 15, nbgn = 15, outgn=1, outsp = 1,sp=1)
ListOut = HGToutCell(ListOut,1)
RES <-Phylter(ListOut, distance="patristic", k=1.5, thres=0.5, quiet=TRUE)

plot.2WR(RES$Complete$mat2WR)

RES$Complete$outgn
RES$Complete$outsp
RES$CellByCell$outcell


##-----------------------------------------------------------------------------------------------------------------
ListOut2 = SimOutliersLg(nbgn =50, nbsp = 50, outgn=4, outsp = 4 ,sp=1)
ListOut2 = BrLengthOutCell(ListOut2,1,9)
RES <-Phylter(ListOut2, distance="patristic", k=1.5, thres=0.5, quiet=TRUE)

plot.2WR(RES$Complete$mat2WR)

RES$Complete$outgn
RES$Complete$outsp
RES$CellByCell$outcell
