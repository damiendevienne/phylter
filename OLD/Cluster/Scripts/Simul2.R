source("/panhome/comte/PhylteR.R")

setwd(dir="/panhome/comte/")

# 2- Analysis
# ***********
tree1 = read.tree("/panhome/comte/ArbreSym.phy")
tree2 = read.tree("/panhome/comte/ArbreAsym.phy")
tree3 = read.tree("/panhome/comte/ArbreRandom.phy")

tr <-c("tree1","tree2","tree3")
nb <- c(100)
outg <- c(0,1,5,10,20,50)
outs <- c(0,1,2,5,15)
outc <- c(0)
outc2 <- c(1,2,5,10,20,50,100)
nrepet = 50

plan1 <- expand.grid(tr=tr, nb = nb, outg = outg, outs = outs, outc = outc)
plan2 <- expand.grid(tr=tr, nb = nb, outg = 0, outs = 0, outc = outc2)
plan=rbind(plan1,plan2)
plan = plan[rep(1:nrow(plan),each = nrepet),]

sp = 1

ListTree = SimOutliersHGT(i_plan=i,tree = get(as.character(plan$tr[i])), nbgn=plan$nb[i], outgn=plan$outg[i], outsp=plan$outs[i], outcell=plan$outc[i], sp=sp)
write.tree(ListTree,paste("/pandata/comte/treeCell/tree",i,sep=""))
