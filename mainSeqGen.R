library(DistatisR)
library(phangorn)
source("/media/aurore/KINGSTON/stageLBBE/pmcoa.R")
source("/home/aurore/Documents/Phylter/Fonctions1.R")
source("/home/aurore/Documents/Phylter/Simulation.R")

setwd(dir="/home/aurore/Documents/Phylter/trees/seqgen/")

tree = rtree(15,rooted = TRUE)
write.tree(tree, file = "arbre.tree")

## Narmol
system("/home/aurore/Téléchargements/Seq-Gen.v1.3.3/source/seq-gen -mHKY -t0.5 -f0.25,0.25,0.25,0.25 -n15 -s0.5 < arbre.tree > seqtrees.dat")
system("phyml -i seqtrees.dat -n 15 -o lr -u arbre.tree --quiet")
ListTrees = read.tree(file="seqtrees.dat_phyml_tree", tree.names = NULL, keep.multi = TRUE)

plot(ListTrees)

##--------------------- Outlier sp-------------------------------------------------------------------------------------------------
##HGT
#1
ListTreesOutSp = HGToutsp(ListTrees,16)
RES <-Phylter(ListTreesOutSp, distance="patristic", k=2, thres=0.5, quiet=TRUE)
RES$Complete$outgn
RES$Complete$outsp
RES$CellByCell$outcell

#2 
tree$edge
ListTreesOutSp2 = list()
for (i in 1:15){
 treeHGT <-HGT(tree,16)
 write.tree(treeHGT, file = "arbreHGT.tree")
  system("/home/aurore/Téléchargements/Seq-Gen.v1.3.3/source/seq-gen -mHKY -t0.5 -f0.25,0.25,0.25,0.25 -n1 -s0.5 < arbreHGT.tree > seqtrees.dat")
  system("phyml -i seqtrees.dat -n 15 -o lr -u arbre.tree --quiet")
  ListTreesOutSp2[[i]]= read.tree(file="seqtrees.dat_phyml_tree")
}
plot(ListTreesOutSp2)

RES <-Phylter(ListTreesOutSp2, distance="patristic", k=2, thres=0.5, quiet=TRUE)
RES$Complete$outgn
RES$Complete$outsp
RES$CellByCell$outcell

##longueur

ListTreesSpLg = BrLengthSp(ListTrees,27)
RES <-Phylter(ListTreesSpLg, distance="patristic", k=2, thres=0.5, quiet=TRUE)
RES$Complete$outgn
RES$Complete$outsp
RES$CellByCell$outcell

plot(ListTreesSpLg)

##--------------------- Outliergn-------------------------------------------------------------------------------------------------