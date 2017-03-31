source("/home/aurore/Documents/Phylter/pmcoa.R")
source("/home/aurore/Documents/Phylter/PhylteR.R")

setwd(dir="/home/aurore/Documents/Phylter/trees/rose/")

ListOut = SimOutliersHGT(nbgn = 10, nbsp = 10, outgn=1, outsp = 1, outcell = 0, sp=1)

RES <-Phylter(ListOut, distance="patristic", k=1.5, thres=0.5, quiet=TRUE, Norm="NONE")
plot.2WR(RES$Complete$mat2WR)

RESCOA <-pMCOA.complete(ListOut, distance="nodal", k=1.5, thres=0.5, quiet=TRUE)
plot.2WR(RESCOA$step1$mat2WR)

RES$Complete$outgn
RES$Complete$outsp
RES$CellByCell$outcell


RESCOA$outcompl$outgn
RESCOA$outcompl$outsp
RESCOA$outcell$outcell

##-----------------------------------------------------------------------------------------------------------------
ListOut2 = SimOutliersLg(nbgn =30, nbsp = 20, outgn=2, outsp = 2, outcell = 1 ,sp=1)

RES <-Phylter(ListOut2, distance="nodal", k=1.5, thres=0.5, quiet=TRUE)
RESCOA <-pMCOA.complete(ListOut2, distance="nodal", k=1.5, thres=0.5, quiet=TRUE)

plot.2WR(RES$Complete$mat2WR)
plot.2WR(RESCOA$step1$mat2WR)

RES$Complete$outgn
RES$Complete$outsp
RES$CellByCell$outcell

RESCOA$outcompl$outgn
RESCOA$outcompl$outsp
RESCOA$outcell$outcell



