##Change un arbre en liste de matrices

trees2matrices.Distatis<-function(trees, distance="nodal",bvalue=0) {
  ##progress bar
  print("Conversion from trees to matrices:")
  pb<-txtProgressBar(style=3, char=".")
  progress<-c(0,1:(length(trees)-1)/(length(trees)-1))
  list.trees<-list()  
  for (i in 1:length(trees)) {
    setTxtProgressBar(pb, progress[i])    
    tree<-trees[[i]]
    if (distance=="nodal") {
      if (bvalue!=0) {
        if (!is.null(tree$node.label)) { 
          l<-1:Nnode(tree)
          indices.nodes<-l[as.numeric(tree$node.label)<bvalue]+Ntip(tree)
          if (length(indices.nodes)>0) {
            for (j in 1:length(indices.nodes)) {
              tree$edge.length[tree$edge[,1]==indices.nodes[j]]<-1e-10
            }   
          }
          tree<-di2multi(tree, tol=1e-9)
        }
        else {
          tree<-di2multi(tree, tol=bvalue)              
        }
      }
      tree.brlen <- compute.brlen(tree, 1)
    }
    else if (distance=="patristic") {
      if (bvalue!=0) {
        l<-1:Nnode(tree)
        indices.nodes<-l[as.numeric(tree$node.label)<bvalue]+Ntip(tree)
        if (length(indices.nodes)>0) {
          for (j in 1:length(indices.nodes)) {
            tree$edge.length[tree$edge[,1]==indices.nodes[j]]<-1e-10
          }   
        }
        tree<-di2multi(tree, tol=1e-9)
      }
      tree.brlen<-tree
    }
    list.trees[[i]]<- tree.brlen
  }
  close(pb)
  TRS<-lapply(list.trees, cophenetic)
  names(TRS)<-as.list(labels(trees))
  return(TRS)
}

#Applique distatis sur une liste de matrice de distance
mat2Dist <-function(matrices){
  genesNumber=length(matrices)
  speciesNumber=nrow(matrices[[genesNumber]])
  TheVeryBigCube = array(0, c(speciesNumber,speciesNumber,genesNumber))
  for (i in 1:genesNumber){
	  TheVeryBigCube[,,i]<-matrices[[i]]
  }
	rownames(TheVeryBigCube) <- rownames(matrices[[i]])
	colnames(TheVeryBigCube) <- colnames(matrices[[i]])
  testDistatis1 <- distatis(TheVeryBigCube)
  dimnames(testDistatis1$res4Splus$PartialF)[[3]] = names(matrices)
  return(testDistatis1)
}

#créé la matrice 2WR a partir des résultats de distatis
Dist2WR <-function(testDistatis2){
  
  Dist.TL = matrix(,dim(testDistatis2$res4Splus$PartialF)[3]*nrow(testDistatis2$res4Splus$PartialF),2)
  colnames(Dist.TL)=c("T","L")
  
  Dist.TLi = matrix(,dim(testDistatis2$res4Splus$PartialF)[3]*nrow(testDistatis2$res4Splus$PartialF),ncol(testDistatis2$res4Splus$PartialF))
  
  colnames(Dist.TLi)=colnames(testDistatis2$res4Splus$PartialF)
  ListRowNames = list()
  j=1
  for (i in 1:dim(testDistatis2$res4Splus$PartialF)[3]) {
    if (j < nrow(Dist.TLi)+1) {
      for (j2 in 1:nrow(testDistatis2$res4Splus$PartialF)){
        a=rownames(testDistatis2$res4Splus$PartialF[,,i])[j2]
        b=dimnames(testDistatis2$res4Splus$PartialF)[[3]][i]
        ListRowNames[[j]] <- paste(a,b,sep=".GEN")
        Dist.TL[j,1]=b
        Dist.TL[j,2]=a
        for (k in 1:ncol(testDistatis2$res4Splus$PartialF[,,i])) {
          Dist.TLi[j,k]=testDistatis2$res4Splus$PartialF[j2,k,i]
        }
        j=j+1
      }
    }
  }
  rownames(Dist.TLi) = ListRowNames
  Dist.TLi=as.data.frame(Dist.TLi)  
  Dist.TL=as.data.frame(Dist.TL) 
  ##Création de la matrice 2WR
  matrixWR<-matrix(ncol=dim(testDistatis2$res4Splus$PartialF)[3], nrow=nrow(testDistatis2$res4Splus$PartialF))
  rownames(matrixWR)=rownames(testDistatis2$res4Splus$PartialF)
  array.vector<-array()
  matrixDif<-matrix(ncol=ncol(Dist.TLi), nrow=nrow(Dist.TLi))
  for (i in 1:nrow(matrixDif)){
    matrixDif[i,]<-(Dist.TLi[i,]-testDistatis2$res4Splus$F[as.character(Dist.TL[i,2]),])^2
    array.vector[i]<-sqrt(sum(matrixDif[i,]))
  }
  colnamesMatrixWR=list()
  for (i in 1:ncol(matrixWR)) {
    matrixWR[,i]<-array.vector[as.character(Dist.TL[,1])==as.character(dimnames(testDistatis2$res4Splus$PartialF)[[3]][i])]
    colnamesMatrixWR[[i]]=dimnames(testDistatis2$res4Splus$PartialF)[[3]][i]
  }
  colnames(matrixWR)=colnamesMatrixWR
  
  return(matrixWR)
}

#Suppression des cell outiers dans les arbres du départ pour ensuite faire la détection des cell outliers

rm.gene.and.species.Distatis<-function(trees, sp2rm, gn2rm) {
  gene.names = list()
  for (i in 1: length(labels(trees))){
    gene.names[i]=labels(trees)[i]
  }
  if (length(CompleteOutlierDist$outsp)>0) {
    for (i in 1:length(trees)) {
      sp<-trees[[i]]$tip.label
      toremove<-intersect(sp, CompleteOutlierDist$outsp)
      trees[[i]]<-drop.tip(trees[[i]], toremove)
    }
  }
  if (length(CompleteOutlierDist$outgn)>0) {
    genes2keep<-setdiff(gene.names, CompleteOutlierDist$outgn)
    trees2<-list()
    j<-0
    for (i in 1:length(trees)) {
      if (is.element(gene.names[i], genes2keep)) {
        j<-j+1
        trees2[[j]]<-trees[[i]]
      }
    }
   names(trees2)<-genes2keep
  }
  else{
    trees2 = trees
  }
  return(trees2)
}

