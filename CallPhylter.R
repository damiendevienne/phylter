require(ape)
require(RSpectra)
require(ggplot2)
require(reshape2)

rm(list=ls())
source("R/trees2matrices.R")
source("R/rename.genes.R")
source("R/impMean.R")
source("R/DistatisFast.R")
source("R/Dist2WR.R")
source("R/normalize.R")
source("R/detect.outliers.R")
source("R/phylter.R")
source("R/simulate.trees.R")
source("R/summary.phylter.R")
source("R/plot.phylter.R")





require(roxygen2)
roxygenize(".",roclets=c("rd","namespace"))








genes<-readLines("AllEukTreeTom")
genenames<-unlist(lapply(strsplit(genes,"/"), function(x) x[length(x)]))
trees<-list()
for (i in 1:length(genes)) {
	trees[[i]]<-read.tree(genes[i])
}
OK<-phylter(trees, gene.names=genenames, thres=0.2)
pdf("dataTom_patristic_thres0.2_k3.pdf",width=20, height=15)
plot(OK, "genes")
plot(OK, "species")
plot2WR(OK)
dev.off()
writeOutput(OK, file="dataTom_patristic_thres0.2_k3.out")
write.tree(fastme.bal(OK$Initial$compromise), file="dataTom_patristic_thres0.2_k3.tre")
write.tree(fastme.bal(OK$Final$compromise), file="dataTom_patristic_thres0.2_k3.tre", append=TRUE)


OK<-phylter(trees, gene.names=genenames, thres=0.2, k=2)
pdf("dataTom_patristic_thres0.2_k2.pdf",width=20, height=15)
plot(OK, "genes")
plot(OK, "species")
plot2WR(OK)
dev.off()
writeOutput(OK, file="dataTom_patristic_thres0.2_k2.out")
write.tree(fastme.bal(OK$Initial$compromise), file="dataTom_patristic_thres0.2_k2.tre")
write.tree(fastme.bal(OK$Final$compromise), file="dataTom_patristic_thres0.2_k2.tre", append=TRUE)




# OK<-phylter(trees, distance="nodal", gene.names=genenames, thres=0.2)
# pdf("dataTom_nodal_thres0.2.pdf",width=20, height=15)
# plot(OK, "genes")
# plot(OK, "species")
# plot2WR(OK)
# dev.off()
# writeOutput(OK, file="dataTom_nodal_thres0.2.out")

# OK<-phylter(trees, distance="nodal", gene.names=genenames, thres=0.2, k=2)
# pdf("dataTom_nodal_thres0.2_k2.pdf",width=20, height=15)
# plot(OK, "genes")
# plot(OK, "species")
# plot2WR(OK)
# dev.off()
# writeOutput(OK, file="dataTom_nodal_thres0.2_k2.out")






ok1<-plot2WR(OK, "Initial")
ok2<-plot2WR(OK2, "Final")



genes<-readLines("ALLNAMES0.5")
trees<-list()
for (i in 1:length(genes)) {
	trees[[i]]<-read.tree(genes[i])
}

# #one dataset
trees<-read.tree("/home/ddevienne/Documents/WORK/Phylter2018/supertree/mamma.tre")
class(trees)<-"multiPhylo"
gene.names<-NULL
trees <- rename.genes(trees, gene.names = gene.names)



trees<-read.tree("ArchaeaTrees.tre")
class(trees)<-"multiPhylo"
gene.names<-NULL
trees <- rename.genes(trees, gene.names = gene.names)



#################
##DONNÉÉS CARINE
#################


folder<-"out_arbre_arbre_online_hmmcleaner_bad_seq_removed"
folder<-"out_arbre_arbre_online_not_clean"
"out_arbre_arbre_online_seuil_0.1"
"out_arbre_arbre_online_seuil_0.2"
"out_arbre_arbre_online_juste_hmmcleaner"
"out_arbre_arbre_online_seuil_0"
"out_arbre_arbre_online_seuil_0.15"
"out_arbre_arbre_online_seuil_0.5"

system(paste("ls datasets/CarineRey/trees/",folder,"/* > CarINEFILES", sep=""))
genes<-readLines("CarINEFILES")
system("rm CarINEFILES")
#read trees
trees<-list()
for (i in 1:length(genes)) {
	trees[[i]]<-read.tree(genes[i])
}
class(trees)<-"multiPhylo"

genenames<-gsub(paste("datasets/CarineRey/trees/",folder,"/", sep=""), "", genes)

OK<-phylter(trees)
