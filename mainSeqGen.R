library(DistatisR)
library(phangorn)

source("/media/aurore/KINGSTON/stageLBBE/pmcoa.R")
source("/home/aurore/Documents/Phylter/Fonctions1.R")
source("/home/aurore/Documents/Phylter/Simulation.R")

setwd(dir="/home/aurore/Documents/Phylter/trees/seqgen/")

tree = rtree(10,rooted = TRUE,min=1,max=10)
#tree = read.tree("arbre.tree")
write.tree(tree, file = "arbre.tree")

## Narmol
system("/home/aurore/Téléchargements/Seq-Gen.v1.3.3/source/seq-gen -mGTR -n20 -s0.5 < arbre.tree > seqtrees.dat")
system("phyml -i seqtrees.dat -n 20 -o lr -u arbre.tree --quiet")
ListTrees = read.tree(file="seqtrees.dat_phyml_tree", tree.names = NULL, keep.multi = TRUE)

RES <-Phylter(ListTrees, distance="patristic", k=2.8, thres=0.5, quiet=TRUE)
RES$Complete$outgn
RES$Complete$outsp
RES$CellByCell$outcell

########################################

ListOut = SimOutliersHGT(nbsp = 5, nbgn = 3, outgn= 0, outsp = 1)
RES <-Phylter(ListOut, distance="patristic", k=2, thres=0.5, quiet=TRUE)
plot.2WR(RES$Complete$mat2WR)
#RES <-pMCOA.complete(ListOut, distance="patristic", k=1.5, thres=0.5, quiet=TRUE)
RES$Complete$outgn
RES$Complete$outsp
RES$CellByCell$outcell

ListOutC = HGToutCell(ListOut, k=1)
RES <-Phylter(ListOutC, distance="patristic", k=2, thres=0.5, quiet=TRUE)
RES$CellByCell$outcell

ListOut2 = SimOutliersLg(nbsp = 20, nbgn = 20, outgn= 2, outsp = 2, sp=1)
RES <-Phylter(ListOut2, distance="patristic", k=2, thres=0.5, quiet=TRUE)
plot.2WR(RES$Complete$mat2WR)
#RES <-pMCOA.complete(ListOut2, distance="patristic", k=1.5, thres=0.5, quiet=TRUE)
RES$Complete$outgn
RES$Complete$outsp
RES$CellByCell$outcell

ListOut2C <- BrLengthOutCell(ListOut2, k=1, ratio=1.5)
RES <-Phylter(ListOut2C, distance="patristic", k=2, thres=0.5, quiet=TRUE)
RES$CellByCell$outcell


