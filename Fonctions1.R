
##----trees2matrices.Distatis change un arbre en liste de matrices----
#Fonction modifiée à partir de celle de pMCOA
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
###On attribut ici de numéro d'un arbre à chaque élément de la liste. (Utile pour la fonction mat2Dist)
  names(TRS)<-as.list(labels(trees))
  return(TRS)
}

##----Gestion des données manquantes dans les matrices----
#Fonction modifiée à partir de celle de pMCOA
gestion.mat.Distatis<-function(matrices) {
  ##progress bar
  print ("Estimation of the quality of the dataset...")
  pb<-txtProgressBar(style=3, char=".")
  progress<-c(0,1:(length(matrices)-1)/(length(matrices)-1))
  qual<-0
  listsp<-colnames(matrices[[1]])
  for (i in 2:length(matrices)) {
    setTxtProgressBar(pb, progress[i])    
    listsp<-union(listsp,colnames(matrices[[i]]))
  }
  length.indiv<-unlist(lapply(lapply(matrices, colnames), length))
  tmp<-rep(1,length(length.indiv))
  testqual<-tmp[(length.indiv==length(listsp))==FALSE]
  if (length(testqual)>0) qual<-1  
  close(pb)
  if (qual==1) {
    pb<-txtProgressBar(style=3, char=".")
    listsp<-colnames(matrices[[1]])
    for (i in 2:length(matrices)) {
      setTxtProgressBar(pb, progress[i])    
      listsp<-union(listsp,colnames(matrices[[i]]))
    }
    close(pb)
    pb<-txtProgressBar(style=3, char=".")
    newcol<-list()
    for (i in 1:length(matrices)) {
      ##print(i)
      setTxtProgressBar(pb, progress[i])    
      newcol[[i]]<-setdiff(listsp,colnames(matrices[[i]]))
      matrices[[i]]<-cbind(matrices[[i]], matrix(ncol=length(newcol[[i]]), nrow=nrow(matrices[[i]]),dimnames=list(colnames(matrices[[i]]),newcol[[i]])))
      matrices[[i]]<-rbind(matrices[[i]], matrix(ncol=ncol(matrices[[i]]), nrow=length(newcol[[i]]),dimnames=list(newcol[[i]],colnames(matrices[[i]]))))
    }
    close(pb)
  }
  mat1 <- matrices[[1]]
  species1 <- row.names(mat1)      
  len2<-length(species1)
  ALL<- list()
  ALL[[1]]<-mat1
  pb<-txtProgressBar(style=3, char=".")  
  for(i in 2:length(matrices)){
    setTxtProgressBar(pb, progress[i])    
    mati <- matrices[[i]]
    speciesi <- row.names(mati)
    indice <- (1:len2)[species1[1] == speciesi]
    for(j in 2:len2){
      indice[j] <- (1:len2)[species1[j] == speciesi]
    }
    mati <- mati[indice,]
    mati <- mati[, indice]
    ALL[[i]] <- mati
  }
  if (qual==1) {
    TEST<-unlist(ALL)
    nbcase<-nrow(ALL[[1]])*nrow(ALL[[1]])
    Nbt<-length(ALL)
    allcase<-0:(Nbt-1)
    MEANS<-array()
    for (i in 1:nbcase) {
      MEANS[i]<-mean(TEST[i+allcase*nbcase],na.rm=TRUE)
    }
    MEANMAT<-matrix(MEANS,nrow=nrow(ALL[[1]]), ncol=ncol(ALL[[1]]))
    rownames(MEANMAT)<-colnames(MEANMAT)<-colnames(ALL[[1]])
    MEANMAT[is.na(MEANMAT)]<-mean(MEANMAT, na.rm=TRUE)
    for (i in 1:length(ALL)) {      
      if (length(newcol[[i]])>0) {
        ALL[[i]][,newcol[[i]]]<-MEANMAT[,newcol[[i]]]
        ALL[[i]][newcol[[i]],]<-MEANMAT[newcol[[i]],]
      }
    }
  }
  close(pb)
  ###On attribut ici de numéro d'un arbre à chaque élément de la liste. (Utile pour la fonction mat2Dist)
  ###A TESTER AVEC UN JEU DE DONNEES A TROUS
  names(ALL)<-names(matrices)
  return(ALL)
}

##----mat2Dist applique distatis sur une liste de matrice de distance----
#En premier lieu, la liste de matrice est changée en cube
mat2Dist <-function(matrices){
  genesNumber=length(matrices)
  speciesNumber=nrow(matrices[[genesNumber]])
  TheVeryBigCube = array(0, c(speciesNumber,speciesNumber,genesNumber))
  for (i in 1:genesNumber){
	  TheVeryBigCube[,,i]<-matrices[[i]]
  }
  #On garde les identifiants de chaque éléments de la matrice dans le cube
	rownames(TheVeryBigCube) <- rownames(matrices[[i]])
	colnames(TheVeryBigCube) <- colnames(matrices[[i]])
	#On applique distatis sur le cube obtenu
  testDistatis1 <- distatis(TheVeryBigCube)
  #On fait suivre le numéro de chaque élement de la matrice de départ aux résultats qui nous interressent
  dimnames(testDistatis1$res4Splus$PartialF)[[3]] = names(matrices)
  return(testDistatis1)
}

###----créé la matrice 2WR à partir des résultats de distatis----
Dist2WR <-function(testDistatis2){
  #Dist.TL = Matrice gène/espèce
  Dist.TL = matrix(,dim(testDistatis2$res4Splus$PartialF)[3]*nrow(testDistatis2$res4Splus$PartialF),2)
  colnames(Dist.TL)=c("T","L")
  #Dist.TLi = Matrice des coordonnées dans les différents axes de tous les couples gène/espèce
  Dist.TLi = matrix(,dim(testDistatis2$res4Splus$PartialF)[3]*nrow(testDistatis2$res4Splus$PartialF),ncol(testDistatis2$res4Splus$PartialF))
  colnames(Dist.TLi)=colnames(testDistatis2$res4Splus$PartialF)
  ListRowNames = list()
  #Remplissage des 2 matrices Dist.TL et Dist.TLi
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
  ##MatrixDif est la matrice des difference au carré entre les coordonnées de chaque couple gène/espèce et les coordonnées moyenne du gène
  matrixDif<-matrix(ncol=ncol(Dist.TLi), nrow=nrow(Dist.TLi))
  for (i in 1:nrow(matrixDif)){
    matrixDif[i,]<-(Dist.TLi[i,]-testDistatis2$res4Splus$F[as.character(Dist.TL[i,2]),])^2
    ##On fait la racine carré de la somme de toutes les différences au carré pour chaque couple gène/espèce et ont les met toutes dans le vecteur suivant
    array.vector[i]<-sqrt(sum(matrixDif[i,]))
  }
  #Remplissage de la matrice 2WR avec les données du vecteurs
  colnamesMatrixWR=list()
  for (i in 1:ncol(matrixWR)) {
    matrixWR[,i]<-array.vector[as.character(Dist.TL[,1])==as.character(dimnames(testDistatis2$res4Splus$PartialF)[[3]][i])]
    #On retrouve le numéro de chaque gène correspondant à chaque distance
    colnamesMatrixWR[[i]]=dimnames(testDistatis2$res4Splus$PartialF)[[3]][i]
  }
  colnames(matrixWR)=colnamesMatrixWR
  return(matrixWR)
}

###-----Suppression des complete outiers dans les arbres du départ pour ensuite faire la détection des cell outliers-----
#Fonction modifiée à partir de celle de pMCOA
rm.gene.and.species.Distatis<-function(trees, sp2rm, gn2rm) {
  gene.names = list()
  for (i in 1: length(labels(trees))){
    gene.names[i]=labels(trees)[i]
  }
  if (length(sp2rm)>0) {
    for (i in 1:length(trees)) {
      sp<-trees[[i]]$tip.label
      toremove<-intersect(sp, sp2rm)
      trees[[i]]<-drop.tip(trees[[i]], toremove)
    }
  }
  if (length(gn2rm)>0) {
    genes2keep<-setdiff(gene.names, gn2rm)
    trees2<-list()
    j<-0
    for (i in 1:length(trees)) {
      if (is.element(gene.names[i], genes2keep)) {
        j<-j+1
        trees2[[j]]<-trees[[i]]
      }
    }
    ##On conserve l'identifiant de chaque gène dans l'arbre de sortie.
   names(trees2)<-genes2keep
  }
  ##Si gn2rm=0, on peut quand même enlever les espèces complete outlier
  else{
    trees2 = trees
  }
  return(trees2)
}

##-------Fonction qui fait tout-----
Fylter <-function(trees, distance="nodal",bvalue=0, k=1, thres=0.5, quiet=TRUE){
  matrices <- trees2matrices.Distatis(trees$trees, distance="nodal",bvalue=0)
  matrices <- gestion.mat.Distatis(matrices)
  Dist <- mat2Dist(matrices)
  WR <- Dist2WR(Dist)
  CompOutl <- detect.complete.outliers(WR, k, thres)
  RES <- NULL
  if (length(CompOutl$outsp)>0 || length(CompOutl$outgn)>0) {
    TREESwithoutCompleteOutlierDist<-rm.gene.and.species.Distatis(trees$trees, CompOutl$outsp, CompOutl$outgn)
    matrices2 = trees2matrices.Distatis(TREESwithoutCompleteOutlierDist, distance ="nodal",bvalue=0)
    matrices2 = gestion.mat.Distatis(matrices2)
    Dist2 <- mat2Dist(matrices2)
    WR2 = Dist2WR(Dist2)
    CellOutl2 <- detect.cell.outliers(WR2, k, quiet)
    RES$Complete <- CompOutl
    RES$CellByCell <- CellOutl2
  }
  else{
    CellOutl2 <- detect.cell.outliers(WR2, k, quiet)
    RES$Complete <- "No Complete Outliers Detected"
    RES$CellByCell <- CellOutl2
  }
  return(RES)
}





