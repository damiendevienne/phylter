library(DistatisR)
library(phangorn)
source("/media/aurore/KINGSTON/stageLBBE/pmcoa.R")
source("/home/aurore/Documents/Phylter/Fonctions1.R")
source("/home/aurore/Documents/Phylter/Simulation.R")


###evolver

setwd(dir="/home/aurore/Documents/Phylter/trees/Evolver/")

tree = rtree(15,rooted = TRUE)
write.tree(tree, file = "arbre.tree")

## Sans outliers
system("perl Write.pl -a arbre.tree -p evolverParam")
system("/home/aurore/Téléchargements/paml4.9c/bin/evolver 5 param.output")
system("phyml -i mc.paml -n 15 -o lr -u arbre.tree --quiet")
ListTreesEvolver = read.tree(file="mc.paml_phyml_tree")

## outlier sp

##HGT
#1
ListTreesEvolverHGT = HGToutsp(ListTreesEvolver,16)
RES <-Phylter(ListTreesEvolverHGT, distance="patristic", k=2, thres=0.5, quiet=TRUE)
RES$Complete$outgn
RES$Complete$outsp
RES$CellByCell$outcell

#2
ListTreesEvolverHGT2 = ListTreesEvolver
for (i in 1:15){
  treeHGT <-HGT(tree,16)
  system("perl Write.pl -a treeHGT -p evolverParamSp")
  system("/home/aurore/Téléchargements/paml4.9c/bin/evolver 5 param.output")
  system("phyml -i mc.paml -n 15 -o lr -u arbre.tree --quiet")
  ListTreesEvolverHGT2[i]= read.tree(file="mc.paml_phyml_tree")
}

##longueur
ListTreesEvolverLG = BrLengthSp(ListTreesEvolver,16)
RES <-Phylter(ListTreesEvolverLG, distance="patristic", k=2, thres=0.5, quiet=TRUE)
RES$Complete$outgn
RES$Complete$outsp
RES$CellByCell$outcell


#write.tree(TREEHGT, file = "arbreHGT.tree")
#system("perl Write.pl -a arbreHGT.tree -p evolverParam")
#system("/home/aurore/Téléchargements/paml4.9c/bin/evolver 5 param.output")
#system("phyml -i mc.paml -n 15 -o lr -u arbre.tree --quiet")
#ListTreesEvolverHGT = read.tree(file="mc.paml_phyml_tree")




ListTrees = ListTrees2

Tree <- ListTrees[[15]]

Treethgt <- HGT(tree,14)
plot(tree)
plot(Treethgt)

treeHGTgn = HGToutgn(tree,15)
plot(tree)
plot(treeHGTgn)

ListTreesHGTsp=HGToutsp(ListTrees2, branche=16)
RES <-Phylter(ListTreesHGTsp, distance="patristic", k=2, thres=0.5, quiet=TRUE)
RES$Complete$outsp

TreeBrLg = BrLength(tree, branche=16, ratio =0.2)
plot(tree)
plot(TreeBrLg)

TreeBrLgGN = BrLengthGn(tree)
plot(tree)
plot(TreeBrLgGN)

TreeBrLgSP = BrLengthSp(ListTrees2, branche=16)
plot(TreeBrLgSP[[5]])
plot(ListTrees2[[5]])

###################################

##Seq-Gen
system("/home/aurore/Téléchargements/Seq-Gen.v1.3.3/source/seq-gen -mHKY -t0.5 -f0.25,0.25,0.25,0.25 -n15 < /home/aurore/Documents/Phylter/trees/arbre.tree > /home/aurore/Documents/Phylter/trees/seqtrees.dat")
system("phyml -i /home/aurore/Documents/Phylter/trees/seqtrees.dat -n 15 -o tlr --quiet")
ListTrees2 = read.tree(file="/home/aurore/Documents/Phylter/trees/seqtrees.dat_phyml_tree", tree.names = NULL, keep.multi = TRUE)

##Rose
ListTreesRose= list()
for (i in 1:15){
  system("rose /home/aurore/Documents/Phylter/trees/roseParam")
  system("phyml -i /home/aurore/Documents/Phylter/trees/RoseTree.phy -n 15 -o tlr --quiet")
  ListTreesRose[[i]] = read.tree(file="/home/aurore/Documents/Phylter/trees/RoseTree.phy_phyml_tree")
}



