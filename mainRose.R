source("/home/aurore/Documents/Phylter/pmcoa.R")
source("/home/aurore/Documents/Phylter/PhylteR.R")

setwd(dir="/home/aurore/Documents/Phylter/trees/rose/")

nbgn = 10
nbsp = 10
outgn=0
outsp = 1
outcell = 0
sp=1

ListOut = SimOutliersHGT(nbgn, nbsp, outgn, outsp, outcell, sp)
write.tree(ListOut,file="test.phy")

mat = trees2matrices.Distatis(ListOut, distance = "patristic")

dist = mat2Dist(mat)

G = dist$res4Cmat$G
rownames(G) = as.character(c(1:nbgn))
GraphDistatisRv(G)

FS=dist$res4Splus$F
PF=dist$res4Splus$PartialF
dimnames(PF)[[3]]=as.character(c(1:nbgn))
GraphDistatisPartial(FS,PF)
GraphDistatisCompromise(FS)


RES <-Phylter(ListOut, distance="patristic", k=1.5, thres=0.5, quiet=TRUE, Norm="NONE", gene.names=as.character(c(1:nbgn)))
plot.2WR(RES$Complete$mat2WR)

RESCOA <-pMCOA.complete(ListOut, distance="nodal", k=1.5, thres=0.5, quiet=TRUE)
plot.2WR(RESCOA$step1$mat2WR)
Dist2WR(mat)

trees2matrices.Distatis(ListOut, distance = "nodal")

RES$Complete$outgn
RES$Complete$outsp
RES$CellByCell$outcell


RESCOA$outcompl$outgn
RESCOA$outcompl$outsp
RESCOA$outcell$outcell

##---------------------------

ListOut = rename.genes(ListOut, gene.names=as.character(c(1:nbgn)))
mat = trees2matrices.Distatis(ListOut, distance = "patristic")
Dist=mat2Dist(mat)
Dist2WR(Dist)


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



