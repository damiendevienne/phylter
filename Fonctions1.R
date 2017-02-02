
##----trees2matrices.Distatis change un arbre en liste de matrices----
trees2matrices.Distatis<-function(trees, distance="nodal", gene.names = NULL) {
  list.trees<-list()  
  for (i in 1:length(trees)) {
    tree<-trees[[i]]
    if (distance=="nodal") {
          tree.brlen <- compute.brlen(tree, 1)
    }
    else if (distance=="patristic") {
          tree.brlen<-tree
    }
    list.trees[[i]]<- tree.brlen
  }
  TRS<-lapply(list.trees, cophenetic)
###On attribut ici de numéro d'un arbre à chaque élément de la liste. (Utile pour la fonction mat2Dist)
  names(TRS)<-as.list(labels(trees))
  return(TRS)
}

#Permet à l'utilisateur d'ajouter des noms de gènes sous forme de liste à la matrice de distance
rename.genes <-function(trees, gene.names=NULL){
  if(is.null(gene.names)==FALSE){
    names(trees)=gene.names
  }
  return(trees)
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
  Distatis <- distatis(TheVeryBigCube)
  #On fait suivre le numéro de chaque élement de la matrice de départ aux résultats qui nous interressent
  dimnames(Distatis$res4Splus$PartialF)[[3]] = names(matrices)
  return(Distatis)
}

##-----Gestion des données manquantes
gestion.matrice<-function(matrices) {
  x = 0
  grandeMatrice=matrix()
  geneNames=list()
  #Création d'une matrice ayant la taille maximale à remplir pour les gènes ayant des espèces manquantes
  for (i in 1: length(matrices)){
   if (nrow(matrices[[i]])>x){
     x=nrow(matrices[[i]])
     grandeMatrice=matrix(nrow=nrow(matrices[[i]]), ncol= ncol(matrices[[i]]), dimnames = dimnames(matrices[[i]]))
   }
  }
  #Création d'une liste contenant les noms des matrices avec des espèces manquantes
  j=1
  for (i in 1: length(matrices)){
    if (nrow(matrices[[i]])<nrow(grandeMatrice)){
      geneNames[[j]]=names(matrices)[[i]]
      j=j+1
    }
  }
  if (length(geneNames)!=0){
  #Pour chaque arbre avec des espèces manquantes, remplissage de la grande matrice avec les données existantes ....
    for (i in 1: length(geneNames)){
      grandeMatrice2 <- grandeMatrice
      for (j in 1: nrow(matrices[[geneNames[[i]]]])){
        for (k in 1: ncol(matrices[[geneNames[[i]]]])){
         grandeMatrice2[j,k]<-matrices[[geneNames[[i]]]][rownames(matrices[[geneNames[[i]]]])[j],colnames(matrices[[geneNames[[i]]]])[k]] 
       }
     }
    # ... puis calcul de la moyenne de tous les arbres complets pour chaque couple gène-espèce manquant
      Inter = setdiff(names(matrices),geneNames)
      for (j in 1: nrow(grandeMatrice2)){
       for (k in 1: ncol(grandeMatrice2)){
          if (is.na(grandeMatrice2[j,k])){
           array = array()
            for (w in 1:length(matrices[Inter])){
              array[[w]] = matrices[Inter][[w]][j,k]
            }
            grandeMatrice2[j,k]=mean(array)
          }
        }
     }
     matrices[[geneNames[[i]]]]=grandeMatrice2
    }
  }
  return(matrices)
}

###----créé la matrice 2WR à partir des résultats de distatis----
Dist2WR <-function(Distatis){
  matrixWR2 = matrix(nrow= length(dimnames(Distatis$res4Splus$PartialF)[[1]]) ,ncol=length(dimnames(Distatis$res4Splus$PartialF)[[3]]))
  colnames(matrixWR2)=dimnames(Distatis$res4Splus$PartialF)[[3]]
  rownames(matrixWR2)=dimnames(Distatis$res4Splus$PartialF)[[1]]
  
  for (i in 1:length(dimnames(Distatis$res4Splus$PartialF)[[3]])){
    for (j in 1: length(dimnames(Distatis$res4Splus$PartialF)[[1]])){
      x = (Distatis$res4Splus$PartialF[dimnames(Distatis$res4Splus$PartialF)[[1]][j],,dimnames(Distatis$res4Splus$PartialF)[[3]][i]]-Distatis$res4Splus$F[dimnames(Distatis$res4Splus$PartialF)[[1]][j],])^2
      matrixWR2[dimnames(Distatis$res4Splus$PartialF)[[1]][j],dimnames(Distatis$res4Splus$PartialF)[[3]][i]] = sqrt(sum(x))
    }
  }
  return(matrixWR2)
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
Phylter <-function(trees, method = "distatis", distance="nodal", k=2, thres=0.5, quiet=TRUE, gene.names=NULL){
  if (is.null(gene.names)==FALSE){
    trees=rename.genes(trees, gene.names=gene.names)
  }
  RES <- NULL
  matrices <- trees2matrices.Distatis(trees, distance=distance)
  #matrices[[1]] = matrices[[1]][1:3,1:3]
  #matrices[[3]] = matrices[[3]][2:5,2:5]
  matrices=gestion.matrice(matrices)
  if (method == "distatis"){
    Dist <- mat2Dist(matrices)
    WR <- Dist2WR(Dist)
    CompOutl <- detect.complete.outliers(WR, k=k, thres=thres)
    if (length(CompOutl$outsp)>0 || length(CompOutl$outgn)>0) {
      TREESwithoutCompleteOutlierDist<-rm.gene.and.species.Distatis(trees, CompOutl$outsp, CompOutl$outgn)
      matrices2 = trees2matrices.Distatis(TREESwithoutCompleteOutlierDist, distance=distance)
      matrices2=gestion.matrice(matrices2)
      Dist2 <- mat2Dist(matrices2)
      WR2 = Dist2WR(Dist2)
      CellOutl2 <- detect.cell.outliers(WR2, k=k, quiet=quiet)
      RES$Complete <- CompOutl
      RES$CellByCell <- CellOutl2
    }
    else{
      CellOutl2 <- detect.cell.outliers(WR,  k=k, quiet=quiet)
      RES$Complete <- CompOutl
      RES$CellByCell <- CellOutl2
    }
  }
  else if (method == "pmcoa"){
    mcoa1 = mat2mcoa(matrices, wtts=NULL, scannf=TRUE, nf="auto")
    WR = mcoa2WRmat(mcoa1)
    colnamesWR = list()
    for (i in 1:length(trees)){
      colnamesWR[[i]]=labels(trees)[[i]]
    }
    colnames(WR)<-colnamesWR
    CompOutl <- detect.complete.outliers(WR, k=k, thres=thres)
    if (length(CompOutl$outsp)>0 || length(CompOutl$outgn)>0) {
      TREESwithoutCompleteOutlierDist<-rm.gene.and.species.Distatis(trees, CompOutl$outsp, CompOutl$outgn)
      newgenenames<-names(TREESwithoutCompleteOutlierDist)
      CellOutl2<-pMCOA(TREESwithoutCompleteOutlierDist, gene.names=newgenenames, distance=distance, bvalue=0, wtts=NULL, scannf=TRUE, nf="auto")
      RES$Complete <- CompOutl
      RES$CellByCell <- CellOutl2
    }
    else{
      CellOutl2 <- detect.cell.outliers(WR, k=k, quiet=quiet)
      RES$Complete <- CompOutl
      RES$CellByCell <- CellOutl2
    }
  }
  return(RES)
}










