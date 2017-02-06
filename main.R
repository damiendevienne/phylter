library(DistatisR)
source("/media/aurore/KINGSTON/stage LBBE/pmcoa.R")
source("/home/aurore/Documents/Phylter/Fonctions1.R")

#création du jeu de données aléatoire: 80 arbres de 60 espèces, avec 15 espèce outlier et 12 genes outliers
trees = gen.trees(Ntrees=20, Ntiptotal=20, Nspmove=3, NbweirdGenes=3)

trees <-read.tree(file="/home/aurore/Téléchargements/Seq-Gen.v1.3.3/seqtrees.dat_phyml_tree")
trees = read.tree(file="/media/aurore/KINGSTON/stage LBBE/Phylter-R/trees/Aguileta-et-al-2008_TREES.txt")

####################################
ListTrees = list()
for (i in 1:100){
  ListTrees[[i]] = trees[[10]]
}
ListTrees[[2]]$edge.length[4]=ListTrees[[2]]$edge.length[4]*2
ListTrees[[45]]$edge.length[9]=ListTrees[[45]]$edge.length[9]/2
ListTrees[[16]]$edge.length[16]=ListTrees[[2]]$edge.length[16]*3
ListTrees[[28]]$edge.length[12]=ListTrees[[45]]$edge.length[12]/2
######################################

RES = AnalyseBranche(trees)
RES$Complete$outgn
RES$Complete$outsp
RES$CellByCell$outcell

#ajout d'outliers
outcell<-add.outliers(trees$trees,3)
trees$trees<-outcell$trees

gene.names=c("gene1", "gene2", "gene3","gene4","gene5","gene6","gene7","gene8","gene9","gene10", "gene11", "gene12", "gene13","gene14","gene15","gene16","gene17","gene18","gene19","gene20")

matrices= trees2matrices.Distatis(trees, distance="nodal")
matrices=gestion.matrice(matrices)
 
#Simulation de données manquantes
Sys1 <- Sys.time()
matrices = gestion.matrice(matrices)
Sys2 <- Sys.time()
Sys1 -Sys2

Sys1 <- Sys.time()
RESdist <-Phylter(trees, distance="patristic", k=1.5, thres=0.4, quiet=TRUE)
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

testDistatis1=mat2Dist(matrices)
Mat2WRDist=Dist2WR(testDistatis1)


## test de Mantel pour voir si les deux matrices 2WR sont corrélées
mantel.test(Mat2WRDist, WRmatMcoa, nrepet = 9999, graph = FALSE)

#recherche des complete outliers pour DISTATIS

CompleteOutlierDist = detect.complete.outliers(matrixWR2,k=2,thres=0.5)
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