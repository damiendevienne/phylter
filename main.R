library(DistatisR)
source("/media/aurore/KINGSTON/stage LBBE/pmcoa.R")
source("/home/aurore/Documents/Phylter/Fonctions1.R")

#création du jeu de données aléatoire: 80 arbres de 60 espèces, avec 15 espèce outlier et 12 genes outliers
trees = gen.trees(Ntrees=50, Ntiptotal=30, Nspmove=8, NbweirdGenes=8)
outcell<-add.outliers(trees$trees, 3)
trees$trees<-outcell$trees
RES <-Fylter(trees, distance="nodal",bvalue=0, k=1, thres=0.5, quiet=TRUE)
RES$Complete$outgn
RES$Complete$outsp
RES$CellByCell$outcell

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

testDistatis1=mat2Dist(matrices)
Mat2WRDist=Dist2WR(testDistatis1)

#recherche des complete outliers pour DISTATIS

CompleteOutlierDist = detect.complete.outliers(Mat2WRDist)
CompleteOutlierDist$outgn
CompleteOutlierDist$outsp

#Recherche des cell outiers pour DISTATIS

TREESwithoutCompleteOutlierDist<-rm.gene.and.species.Distatis(trees$trees, CompleteOutlierDist$outsp, CompleteOutlierDist$outgn)

matrices2 = trees2matrices.Distatis(TREESwithoutCompleteOutlierDist, distance ="nodal",bvalue=0)
matrices2 = gestion.mat.Distatis(matrices2)
testDistatis2 <- mat2Dist(matrices2)
Mat2WRDist2 = Dist2WR(testDistatis2)
plot.2WR(Mat2WRDist)
plot.2WR(Mat2WRDist2)

OutlerCell <- detect.cell.outliers(Mat2WRDist2)