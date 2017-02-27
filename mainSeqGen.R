library(DistatisR)
library(phangorn)
source("/media/aurore/KINGSTON/stageLBBE/pmcoa.R")
source("/home/aurore/Documents/Phylter/Fonctions1.R")
source("/home/aurore/Documents/Phylter/Simulation.R")

setwd(dir="/home/aurore/Documents/Phylter/trees/seqgen/")

tree = rtree(30,rooted = TRUE)
#tree = read.tree("arbre.tree")
write.tree(tree, file = "arbre.tree")

## Narmol
system("/home/aurore/Téléchargements/Seq-Gen.v1.3.3/source/seq-gen -mGTR -n20 -s0.5 < arbre.tree > seqtrees.dat")
system("phyml -i seqtrees.dat -n 20 -o lr -u arbre.tree --quiet")
ListTrees = read.tree(file="seqtrees.dat_phyml_tree", tree.names = NULL, keep.multi = TRUE)

plot(ListTrees)

RES <-Phylter(ListTrees, distance="patristic", k=2.8, thres=0.5, quiet=TRUE)
RES$Complete$outgn
RES$Complete$outsp
RES$CellByCell$outcell

##--------------------- Outlier sp-------------------------------------------------------------------------------------------------
##HGT
#1
ListTreesOutSp = HGToutsp(ListTrees,14)
RES <-Phylter(ListTreesOutSp, distance="patristic", k=2.8, thres=0.5, quiet=TRUE)
RES$Complete$outgn
RES$Complete$outsp
RES$CellByCell$outcell

ListTrees[[1]]$edge

#2 
tree$edge
ListTreesOutSp2 = list()
for (i in 1:20){
 treeHGT <-HGT(tree,14)
 write.tree(treeHGT, file = "arbreHGT.tree")
  system("/home/aurore/Téléchargements/Seq-Gen.v1.3.3/source/seq-gen -mGTR -n1 -d1 < arbre.tree > seqtrees.dat")
  system("phyml -i seqtrees.dat -n 1 -o lr -u arbre.tree --quiet")
  ListTreesOutSp2[[i]]= read.tree(file="seqtrees.dat_phyml_tree")
}
for (i in 1:length(ListTreesOutSp2)){
  plot(ListTreesOutSp2[[i]])
}

RES <-Phylter(ListTreesOutSp2, distance="patristic", k=2.8, thres=0.5, quiet=TRUE)
RES$Complete$outgn
RES$Complete$outsp
RES$CellByCell$outcell

##longueur

#1
ListTreesSpLg = BrLengthSp(ListTrees,16)
RES <-Phylter(ListTreesSpLg, distance="patristic", k=2, thres=0.5, quiet=TRUE)
RES$Complete$outgn
RES$Complete$outsp
RES$CellByCell$outcell

#2
#2 
ListTreesOutSpLg2 = list()
for (i in 1:20){
  treeLG <-BrLength(tree,14, ratio=5)
  write.tree(treeLG, file = "arbreHGT.tree")
  system("/home/aurore/Téléchargements/Seq-Gen.v1.3.3/source/seq-gen -mGTR -n1 -s0.5 < arbre.tree > seqtrees.dat")
  system("phyml -i seqtrees.dat -n 1 -o lr -u arbre.tree --quiet")
  ListTreesOutSpLg2[[i]]= read.tree(file="seqtrees.dat_phyml_tree")
}
RES <-Phylter(ListTreesOutSpLg2, distance="patristic", k=2.8, thres=0.5, quiet=TRUE)
RES$Complete$outgn
RES$Complete$outsp
RES$CellByCell$outcell

##--------------------- Outliergn-------------------------------------------------------------------------------------------------
##HGT

ListTreesGpHGT = HGToutgn(ListTrees, k=3)
RES <-Phylter(ListTreesGpHGT, distance="patristic", k=2, thres=0.5, quiet=TRUE)
RES$Complete$outgn
RES$Complete$outsp
RES$CellByCell$outcell


##longueur

ListTreesLgn = BrLengthGn(ListTrees, Listgn=c(3,6,7), ratioMin = 5, ratioMax = 10)
RES <-Phylter(ListTreesLgn, distance="patristic", k=2, thres=0.5, quiet=TRUE)
RES$Complete$outgn
RES$Complete$outsp
RES$CellByCell$outcell

##--------------------- Outliercell-------------------------------------------------------------------------------------------------
##HGT



##longueur

ListTreesLoutcell = BrLengthOutCell(ListTrees, n=1, ratio=3)
RES <-Phylter(ListTreesLoutcell, distance="patristic", k=2, thres=0.5, quiet=TRUE)
RES$Complete$outgn
RES$Complete$outsp
RES$CellByCell$outcell

########################################""""

ListOut = SimOutliersHGT(10,10,1,1)
RES <-Phylter(ListOut, distance="patristic", k=2, thres=0.5, quiet=TRUE)
RES$Complete$outgn
RES$Complete$outsp
RES$CellByCell$outcell

ListOutC = HGToutCell(ListOut, n=1)
RES <-Phylter(ListOutC, distance="patristic", k=2, thres=0.5, quiet=TRUE)
RES$CellByCell$outcell

ListOut2 = SimOutliersLg(20,20,2,2)
RES <-Phylter(ListOut2, distance="patristic", k=2, thres=0.5, quiet=TRUE)
RES$Complete$outgn
RES$Complete$outsp
RES$CellByCell$outcell

ListOut2C <- BrLengthOutCell(ListOut2, n=1)
RES <-Phylter(ListOut2C, distance="patristic", k=2, thres=0.5, quiet=TRUE)
RES$CellByCell$outcell
