library(DistatisR)
library(phangorn)
source("/media/aurore/KINGSTON/stageLBBE/pmcoa.R")
source("/home/aurore/Documents/Phylter/Fonctions1.R")
source("/home/aurore/Documents/Phylter/Simulation.R")

tree = rtree(15,rooted = TRUE)
write.tree(tree, file = "/home/aurore/Documents/Phylter/trees/arbre.tree")

system("/home/aurore/Téléchargements/Seq-Gen.v1.3.3/source/seq-gen -mHKY -t0.5 -f0.25,0.25,0.25,0.25 -n15 < /home/aurore/Documents/Phylter/trees/arbre.tree > /home/aurore/Documents/Phylter/trees/seqtrees.dat")

system("phyml -i /home/aurore/Documents/Phylter/trees/seqtrees.dat -n 15 -o tlr --quiet")

ListTrees2 = read.tree(file="/home/aurore/Documents/Phylter/trees/seqtrees.dat_phyml_tree", tree.names = NULL, keep.multi = TRUE)

RES <-Phylter(ListTrees, distance="patristic", k=2, thres=0.5, quiet=TRUE)
RES$Complete$outgn
RES$Complete$outsp

ListTrees = ListTrees2

Tree <- ListTrees[[15]]

Treethgt <- HGT(Tree,14)
plot(Tree)
plot(Treethgt)

