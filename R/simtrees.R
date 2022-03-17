# Simple(simplistic) simulation of trees with outliers. No missing data. 

#' Simplistic simulation of gene trees with outliers
#' 
#' Simple (simplistic) simulation of trees with outliers. 
#' 
#' The simulation process is as follows: a first tree is generated with the \code{rtree()} function
#' and is then duplicated and modified according to the parameters chosen by the user.
#' 
#' @usage simtrees(Ngn, Nsp, Nsp.out = 0, Ngn.out = 0, Nb.cell.outlier = 0, brlen.sd = 0)
#' @param Ngn Number of gene trees to simulate. 
#' @param Nsp Number of species (tips) per tree.
#' @param Nsp.out Number of outlier species (also called rogue taxa). 0 = none.
#' @param Ngn.out Number of outlier genes. 0 = none.
#' @param Nb.cell.outlier Number of times one species in one tree is misplaced. 0 = none.
#' @param brlen.sd Heterogeneity of branch lengths in trees. A value with mean 0 and standard
#' deviation equal to brlen.sd is added to each branch length.
#' @return A list of trees in \code{multiPhylo} format.
#' @examples  
#' # Very basic simulator, for debuggin purpose mainly.
#' # examples: 30 genes, 120 species, 2 outlier species, 3 outlier genes
#' # 4 gene/species outliers, branch length variance = 0.6 
#' trees<-simtrees(30,120,2,3,4,0.6)
#' 
#' 
#' 
#' @importFrom ape rtree drop.tip Ntip bind.tree multi2di
#' @importFrom stats rnorm
#' @export
simtrees<-function(Ngn, Nsp, Nsp.out=0,Ngn.out=0,Nb.cell.outlier=0, brlen.sd=0) {
  gen.trees<-function(Ntrees, Ntiptotal, Nspmove, NbweirdGenes, brlen.sd) {
    res<-list()
    selec<-NULL
    anames<-NULL
    a<-rtree(Ntiptotal)
    TrueTree<-a
    if (Nspmove>0) {
      anames<-sample(a$tip.label, Nspmove)
      a<-drop.tip(a, anames)
    }
    for (i in 1:Ntrees){
      res[[i]]<-a
      if (Nspmove>0) {
        for (j in 1:Nspmove) {
          b<-rtree(2)
          nodelabs<-(Ntip(res[[i]])+1):(Ntip(res[[i]])+Ntip(res[[i]])-1)
          tiplab<-c(anames[j], "out")
          b$tip.label<-tiplab
          tr<-bind.tree(res[[i]],b, where=sample(nodelabs, 1))
          tr<-drop.tip(tr, tiplab[2])
          tr<-multi2di(tr)
          tr$edge.length[tr$edge.length==0]<-mean(tr$edge.length)
          res[[i]]<-tr
        }
      }
    }
    if (NbweirdGenes>0) {
      sam<-1:Ntrees
      selec<-sample(sam,NbweirdGenes)
      for (k in selec) {
        res[[k]]<-rtree(Ntiptotal)
      }
    }
    RES<-NULL
    RES$TrueTree<-TrueTree
    RES$trees<-lapply(res, SlightChangeInTreeBrLen, brlen.sd=brlen.sd)
    RES$gn<-selec  
    RES$sp<-anames
    return(RES)
  }

  add.outliers<-function(trees, nbrep) {
    move.tip<-function(tree, tip) {
      treesmall<-drop.tip(tree, tip)
      temptree<-rtree(2)
      temptree$tip.label[1]<-"out"
      temptree$tip.label[2]<-tip
      node.attach<-sample(unique(treesmall$edge[,1]),1)
      treenew<-bind.tree(treesmall, temptree,where=node.attach)
      treenew<-drop.tip(treenew,"out")
      treenew<-multi2di(treenew)
      return(treenew)
    }
    Ntiptotal<-Ntip(trees[[1]])
    Ntrees<-length(trees)
    rrr<-c(NA,NA)
    names(rrr)<-c("Species","Genes")
    for (w in 1:nbrep) {
      which.gn<-sample(1:Ntrees, 1)
      which.sp<-paste("t",sample(1:Ntiptotal, 1),sep="")
      rrr<-rbind(rrr,c(which.sp,which.gn))
      trees[[which.gn]]<-move.tip(trees[[which.gn]], which.sp)
    }
    RES<-NULL
    RES$trees<-trees
    RES$outl<-rrr[2:nrow(rrr),]  
    return(RES)
  }
  SlightChangeInTreeBrLen<-function(tr, brlen.sd) {
    tr$edge.length<-abs(tr$edge.length+rnorm(length(tr$edge.length), sd=brlen.sd))
    tr
  }
  tr<-gen.trees(Ngn,Nsp,Nsp.out,Ngn.out, brlen.sd)
  if (Nb.cell.outlier>0) trees<-add.outliers(tr$trees, Nb.cell.outlier)$trees
  else trees<-tr$trees
  return(trees)
}
