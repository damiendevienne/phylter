library(DistatisR)
source("/media/aurore/KINGSTON/stage LBBE/pmcoa.R")
source("/home/aurore/Documents/Phylter/Fonctions1.R")

#création du jeu de données aléatoire: 80 arbres de 60 espèces, avec 15 espèce outlier et 12 genes outliers
trees = gen.trees(Ntrees=200, Ntiptotal=150, Nspmove=9, NbweirdGenes=10)

#ajout d'outliers
outcell<-add.outliers(trees$trees,5)
trees$trees<-outcell$trees

matrices = trees2matrices.Distatis(trees$trees)

Sys1 <- Sys.time()
RESdist <-Phylter(trees$trees, method="distatis", distance="nodal", k=1.2, thres=0.4, quiet=TRUE)
Sys2 <- Sys.time()
Sys1 -Sys2

Sys1 <- Sys.time()
RESpmcoa <-Phylter(trees$trees, method="pmcoa", distance="nodal", k=1.2, thres=0.4, quiet=TRUE)
Sys2 <- Sys.time()
Sys1 -Sys2

RESdist$Complete$outgn
RESdist$Complete$outsp
RESdist$CellByCell$outcell

#-------------------------------------------------------------------------------------------------------
#PHYLOMCOA

mcoa=mat2mcoa(matrices)
WRmatMcoa = mcoa2WRmat(mcoa)
#CompleteOutlierMcoa = detect.complete.outliers(WRmatMcoa)
plot.2WR(WRmatMcoa)
pMCOAout=pMCOA.complete(trees$trees)

#recherche des complete outliers pour MCOA

pMCOAout$outcompl$outgn
pMCOAout$outcompl$outsp

#Recherche des cell outiers pour MCOA

pMCOAout$outcell$outcell

#----------------------------------------------------------------------------------------------
#DISTATIS

testDistatis2=mat2Dist(matrices)
Mat2WRDist=Dist2WR(testDistatis1)

#recherche des complete outliers pour DISTATIS

CompleteOutlierDist = detect.complete.outliers(matrixWR2,k=0.8,thres=0.2)
CompleteOutlierDist$outgn
CompleteOutlierDist$outsp

#Recherche des cell outiers pour DISTATIS

TREESwithoutCompleteOutlierDist<-rm.gene.and.species.Distatis(trees$trees, CompleteOutlierDist$outsp, CompleteOutlierDist$outgn)

matrices2 = trees2matrices.Distatis(TREESwithoutCompleteOutlierDist, distance ="nodal")
testDistatis2 <- mat2Dist(matrices2)
Mat2WRDist2 = Dist2WR(testDistatis2)
plot.2WR(Mat2WRDist)
plot.2WR(Mat2WRDist2)

OutlerCell <- detect.cell.outliers(Mat2WRDist2)