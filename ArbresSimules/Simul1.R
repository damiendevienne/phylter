source("/home/aurore/Documents/Phylter/pmcoa.R")
source("/home/aurore/Documents/Phylter/PhylteR.R")

setwd(dir="/home/aurore/Documents/Phylter/ArbresSimules/")

#Arguments
args <- commandArgs(TRUE)
tree <- as.character(args[1])
nbgn <- as.integer(args[2])
outgn <- as.integer(args[3])
outsp <- as.integer(args[4])
outcell <- as.integer(args[5])
i <- as.character(args[6])

#parametres immobiles
tr = read.tree(tree)
sp = 1

system(paste("echo", i, ">> log", sep= " "))

ListTree = SimOutliersHGT(tr, nbgn, outgn, outsp, outcell, sp)

write.tree(ListTree,paste("/home/aurore/Documents/Phylter/ArbresSimules/arbre/arbre",i))

system("echo *** >> log")
