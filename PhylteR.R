##load required packages
require(missMDA)
require(FactoMineR)
require(DistatisR)
require(ape)
require(phangorn)
require(adephylo)

##----trees2matrices.Distatis change un arbre en liste de matrices----
trees2matrices.Distatis<-function(trees, distance="nodal") {
  list.trees<-list()  
  for (i in 1:length(trees)) {
    tree<-trees[[i]]
    if (distance=="nodal") {
      tree.brlen <- distTips(trees[[i]], method="nNodes")
    }
    else if (distance=="patristic") {
      tree.brlen<- distTips(trees[[i]],method="patristic")
    }
    list.trees[[i]]<-as.matrix(tree.brlen)
  }
  if (!is.null(names(trees))){
    names(list.trees) = names(trees)
  }
  return(list.trees)
}

#Permet à l'utilisateur d'ajouter des noms de gènes (sous forme de liste) aux arbres
rename.genes <-function(trees, gene.names=NULL){
  if(!is.null(gene.names)){
    names(trees)=gene.names
  }
  else{
    if(is.null(names(trees))){
      names(trees) = as.character(c(1:length(trees)))
    }
  }
  return(trees)
}

##----mat2Dist applique distatis sur une liste de matrices de distance----
##Norm = "none" si on ne veut pas normaliser les données --> permet de detecter les variations de longueur dans des gènes entiers
##Norm = "mfa" normalise les données
mat2Dist <-function(matrices, Norm="NONE"){
  #En premier lieu, la liste de matrice est changée en cube
  row = rownames(matrices[[1]])
  for (i in 1:length(matrices)){
    matrices[[i]]=matrices[[i]][row,row]
  }
  genesNumber=length(matrices)
  speciesNumber=nrow(matrices[[genesNumber]])
  TheVeryBigCube = array(0, c(speciesNumber,speciesNumber,genesNumber))
  for (i in 1:genesNumber){
    TheVeryBigCube[,,i]<-matrices[[i]]
  }
  rownames(TheVeryBigCube) <- rownames(matrices[[1]])
  colnames(TheVeryBigCube) <- colnames(matrices[[1]])
  #On applique distatis sur le cube obtenu
  Distatis <- distatis(TheVeryBigCube, Norm = Norm)
  #On fait suivre le numéro de chaque élement de la matrice de départ aux résultats qui nous interressent
  dimnames(Distatis$res4Splus$PartialF)[[3]] = names(matrices)
  return(Distatis)
}

##Fonction utilisant le package missMDA afin dimputer des données d'espèces potentiellement manquantes dans les différentes matrices
impPCA.multi<-function(matrices, ncp = 3, center =FALSE, scale = FALSE, maxiter = 10000) {
  species = list()
  grandeMatrice=matrix()
  geneNames=list()
  #Création d'une matrice ayant la taille maximale à remplir pour les gènes ayant des espèces manquantes
  for (i in 1: length(matrices)){
    species = c(species,setdiff(rownames(matrices[[i]]),species))
  }
  grandeMatrice=matrix(nrow=length(species), ncol= length(species))
  rownames(grandeMatrice) = species
  colnames(grandeMatrice) = species
  #Création d'une liste contenant les noms des matrices avec des espèces manquantes
  j=1
  for (i in 1: length(matrices)){
    if (nrow(matrices[[i]])<nrow(grandeMatrice)){
      if(!is.null(names(matrices))){
        geneNames[j]=names(matrices)[[i]]
      }
      else{
        geneNames[j]=i
      }
      j=j+1
    }
  }
  if (length(geneNames)!=0){
    #Pour chaque arbre avec des espèces manquantes, remplissage de la grande matrice avec les données existantes ....
    matrices2=matrices
    for (i in 1: length(geneNames)){
      row = rownames(matrices2[[geneNames[[i]]]])
      grandeMatrice2 <- grandeMatrice
      grandeMatrice2[row,row]<-matrices2[[geneNames[[i]]]][row,row]
      matrices2[[geneNames[[i]]]]<-grandeMatrice2
    }
    matrices3=matrices2
    row = rownames(matrices2[[1]])
    for (i in 1:length(matrices2)){
      matrices3[[i]]=as.vector(matrices2[[i]][row,row][upper.tri(matrices2[[i]][row,row])])
    }
    mat = do.call(cbind,matrices3)
    matIPCA = imputePCA(mat, center =center, scale = scale, maxiter = maxiter) ##estimation des données manquantes par ACP itérative
    matIPCA <- matIPCA$completeObs
    matricesFT=list()
    for (i in 1:length(matrices2)){
      matrices[[i]]= matrix(nrow=length(row),ncol=length(row))
      rownames(matrices[[i]])=row
      colnames(matrices[[i]])=row
      data=as.vector(matIPCA[,i])
      matrices[[i]][upper.tri(matrices2[[i]])]=data
      matricesFT[[i]]=t(matrices[[i]])
      matrices[[i]][lower.tri(matrices2[[i]])]=matricesFT[[i]][lower.tri(matrices2[[i]])]
      for (j in 1:ncol(matrices[[i]])){
        matrices[[i]][j,j] = 0
      }
    }
  }
  return(matrices)
}

#############IMPUTE PCA#################
imputePCA <- function (X, ncp = 2, center= FALSE, scale=FALSE, method=c("Regularized","EM"),row.w=NULL,coeff.ridge=1,threshold = 1e-6,seed = NULL,nb.init=1,maxiter=1000,...){
  impute <- function (X, ncp = 4, center=FALSE, scale=FALSE, method=NULL,threshold = 1e-6,seed = NULL,init=1,maxiter=1000,row.w=NULL,coeff.ridge=1,...){
    moy.p <- function(V, poids) {
      res <- sum(V * poids,na.rm=TRUE)/sum(poids[!is.na(V)])
    }
    ec <- function(V, poids) {
      res <- sqrt(sum(V^2 * poids,na.rm=TRUE)/sum(poids[!is.na(V)]))
    }
    nb.iter <- 1
    old <- Inf
    objective <- 0
    if (!is.null(seed)){set.seed(seed)}
    X <- as.matrix(X)
    ncp <- min(ncp,ncol(X),nrow(X)-1)
    missing <- which(is.na(X))
    mean.p <- apply(X, 2, moy.p,row.w)
    if (center!= TRUE){
      mean.p=rep(0,ncol(X))
    }
    Xhat <- t(t(X)-mean.p)
    et <- apply(Xhat, 2, ec,row.w)
    if (scale) Xhat <- t(t(Xhat)/et)
    if (any(is.na(X))) Xhat[missing] <- 0
    if (init>1) Xhat[missing] <- rnorm(length(missing)) ## random initialization
    fittedX <- Xhat
    if (ncp==0) nb.iter=0
    while (nb.iter > 0) {
      Xhat[missing] <- fittedX[missing]
      if (scale) Xhat=t(t(Xhat)*et)
      Xhat <- t(t(Xhat)+mean.p)
      if (center== TRUE){
        mean.p <- apply(Xhat, 2, moy.p,row.w)
      }
      Xhat <- t(t(Xhat)-mean.p)
      et <- apply(Xhat, 2, ec,row.w)
      if (scale) Xhat <- t(t(Xhat)/et)
      svd.res <- svd.triplet(Xhat,row.w=row.w,ncp=ncp)
      #       sigma2 <- mean(svd.res$vs[-(1:ncp)]^2)
      sigma2  <- nrow(X)*ncol(X)/min(ncol(X),nrow(X)-1)* sum((svd.res$vs[-c(1:ncp)]^2)/((nrow(X)-1) * ncol(X) - (nrow(X)-1) * ncp - ncol(X) * ncp + ncp^2))
      sigma2 <- min(sigma2*coeff.ridge,svd.res$vs[ncp+1]^2)
      if (method=="em") sigma2 <-0
      lambda.shrinked=(svd.res$vs[1:ncp]^2-sigma2)/svd.res$vs[1:ncp]
      fittedX = tcrossprod(t(t(svd.res$U[,1:ncp,drop=FALSE]*row.w)*lambda.shrinked),svd.res$V[,1:ncp,drop=FALSE])
      fittedX <- fittedX/row.w
      ## Pas de valeurs négatives
      fittedX[fittedX<0]=0
      diff <- Xhat-fittedX
      diff[missing] <- 0
      objective <- sum(diff^2*row.w)
      #       objective <- mean((Xhat[-missing]-fittedX[-missing])^2)
      criterion <- abs(1 - objective/old)
      old <- objective
      nb.iter <- nb.iter + 1
      if (!is.nan(criterion)) {
        if ((criterion < threshold) && (nb.iter > 5))  nb.iter <- 0
        if ((objective < threshold) && (nb.iter > 5))  nb.iter <- 0
      }
      if (nb.iter > maxiter) {
        nb.iter <- 0
        warning(paste("Stopped after ",maxiter," iterations"))
      }
    }
    if (scale) Xhat <- t(t(Xhat)*et)
    Xhat <- t(t(Xhat)+mean.p)
    completeObs <- X
    completeObs[missing] <- Xhat[missing]
    if (scale) fittedX <- t(t(fittedX)*et)
    fittedX <- t(t(fittedX)+mean.p)
    
    result <- list()
    result$completeObs <- completeObs
    result$fittedX <- fittedX
    return(result) 
  }
  
  #### Main program
  method <- match.arg(method,c("Regularized","regularized","EM","em"),several.ok=T)[1]
  obj=Inf
  method <- tolower(method)
  if (ncp>min(nrow(X)-2,ncol(X)-1)) stop("ncp is too large")
  if (is.null(row.w)) row.w = rep(1,nrow(X))/nrow(X)
  for (i in 1:nb.init){
    if (!any(is.na(X))) return(X)
    res.impute=impute(X, ncp=ncp, scale=scale,center= center, method=method, threshold = threshold,seed=if(!is.null(seed)){(seed*(i-1))}else{NULL},init=i,maxiter=maxiter,row.w=row.w,coeff.ridge=coeff.ridge)
    if (mean((res.impute$fittedX[!is.na(X)]-X[!is.na(X)])^2) < obj){
      res <- res.impute
      obj <- mean((res.impute$fittedX[!is.na(X)]-X[!is.na(X)])^2)
    }
  }
  return(res)
}


###----créé la matrice 2WR à partir des résultats de distatis----
Dist2WR <-function(Distatis){
  matrixWR2 = matrix(nrow= dim(Distatis$res4Splus$PartialF)[[1]] ,ncol=dim(Distatis$res4Splus$PartialF)[[3]])
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
Phylter <-function(trees, distance="nodal", k=2, thres=0.5, quiet=TRUE, gene.names=NULL, Norm="NONE"){
  trees=rename.genes(trees, gene.names=gene.names)
  RES <- NULL
  matrices <- trees2matrices.Distatis(trees, distance=distance)
  matrices=impPCA.multi(matrices)
  Dist <- mat2Dist(matrices,Norm=Norm)
  WR <- Dist2WR(Dist)
  CompOutl <- detect.complete.outliers(WR, k=k, thres=thres)
  if (length(CompOutl$outsp)>0 || length(CompOutl$outgn)>0) {
    TREESwithoutCompleteOutlierDist<-rm.gene.and.species.Distatis(trees, CompOutl$outsp, CompOutl$outgn)
    matrices2 = trees2matrices.Distatis(TREESwithoutCompleteOutlierDist, distance=distance)
    matrices2=impPCA.multi(matrices2)
    Dist2 <- mat2Dist(matrices2,Norm=Norm)
    WR2 = Dist2WR(Dist2)
    CellOutl2 <- detect.cell.outliers(WR2, k=k+2, quiet=quiet)
    RES$Complete <- CompOutl
    RES$CellByCell <- CellOutl2
  }
  else{
    CellOutl2 <- detect.cell.outliers(WR,  k=k+2, quiet=quiet)
    RES$Complete <- CompOutl
    RES$CellByCell <- CellOutl2
  }
  return(RES)
}

######-----------------------------------Simulation--------------------------------------------------------------------------

###Fonction qui génère une liste d'arbres de nbgn gènes avec nbsp espèces contenant des outliers gènes (outgn) et espèces (outsp) générés par HGT
##nbsp = nombre d'espèces dans l'arbre / nbgn = nombre d'arbres / outgn = nb d'outlier gènes /outsp = nb d'oulier sp  /outcell = couple outlier gn/sp
##sp= 0 -> HGT sur tous les branches / sp = 1 -> HGT que sur les branches extérieures
SimOutliersHGT <-function(tree, nbgn, nbsp, outgn, outsp, outcell, sp = 0){
  genes=NULL
  #tree<-rtree(nbsp,rooted = TRUE,min=1,max=10)
  #write.tree(tree, file = "arbre.tree")
  ListOutGnTree =list()
  "multiPhylo"->class(ListOutGnTree)
  for (i in 1:nbgn){ # n fois le même arbre qui va evoluer différement
    ListOutGnTree[[i]]=tree
  }
  if(outsp!=0){
    if(sp==1){
      ## On ne fait des hgt que sur les branches externes (pour pouvoir contrôler le nombre d'espèces outliers)
      samp = sample(tree$tip.label,outsp) 
      print(paste("outsp",samp,sep=" "))
      species = samp
      for (o in 1:outsp){
        ListOutGnTree = HGToutsp(ListOutGnTree, samp[o])
      }
    }
    else if(sp==0){
      #C'est possible que deux branches différentes choisies comme outliers aient les mêmes espèces déscendantes
      bra = c(1:nrow(tree$edge))
      br<-sample(bra,outsp)
      for (B in 1:length(br)){
        bran<-tree$edge[br[B],2]
        samp<-tree$tip.label[Descendants(tree, bran, "tips")[[1]]]
        print(paste("outsp",samp,sep=" "))
        species = samp
        ListOutGnTree = HGToutsp(ListOutGnTree,samp)
      }
    }
  }
  else {species = NULL}
  if (outgn !=0){
    samp = sample(1:length(ListOutGnTree),outgn)
    print(paste("outgn",samp,sep=" "))
    genes = samp
    for (j in 1:length(samp)){
      ListOutGnTree[[samp[j]]] =  HGT1gn(ListOutGnTree[[samp[j]]], n=length(ListOutGnTree[[samp[j]]]$edge))
    }
  }
  else{genes = NULL}
  if (outcell !=0){ ## Les outcell ne peuvent pas être sur les espèces ou sur les gènes déjà complete outlier. 
    #Si on veut que ce soit possible d'avoir des outcell sur des gènes ou sp outliers complets, on fixe outcell = 0 et on créé des outcells avec la fonction
    #HGToutCell()
    if (!is.null(genes)){genespossible = c(1:length(ListOutGnTree))[which(c(1:length(ListOutGnTree))!=genes)]}
    else{genespossible=c(1:length(ListOutGnTree))}
    samp= sample(genespossible,outcell)
    for (i in 1:length(samp)){
      if (!is.null(species)){speciespossible = setdiff(ListOutGnTree[[samp[i]]]$tip.label, species)}
      else{speciespossible=ListOutGnTree[[samp[i]]]$tip.label}
      lab = sample(speciespossible,1)
      print(paste("outcell =", samp[i],"/", lab[[1]],sep=" "))
      ListOutGnTree[[samp[i]]] = HGT2(ListOutGnTree[[samp[i]]],lab[[1]])
    }
  }
  ListTreesOut=list()
  "multiPhylo"->class(ListTreesOut)
  compteur= 0
  for  (i in 1: length(ListOutGnTree)){
    tree<-di2multi(ListOutGnTree[[i]], tol=1e-3)
    n=length(tree$tip.label)
    write.tree(tree, file = "arbreHGT.tree")
    system(paste("perl WriteRose.pl -a arbreHGT.tree -p roseParam -l ", n, sep=""))
    system("rose param.output")
    system("phyml -i RoseTree.phy -m JC69 -n 1 -o tlr --quiet > sortie.out")
    ListTreesOut[[i]]= read.tree(file="RoseTree.phy_phyml_tree")
    compteur=compteur+1
    print(compteur)
  }
  #RES=NULL
  #RES$genes=genes
  #RES$ListTreesOut=ListTreesOut
  #return(RES)
  #return(ListOutGnTree)
  return(ListTreesOut)
}

###Fonction qui génère une liste d'arbres de nbgn gènes avec nbsp espèces contenant des outliers gènes (outgn) et espèces (outsp) générés
###en modifiant les longueurs de branches
#nbsp = nombre d'espèces dans l'arbre / nbgn = nombre d'arbres / outgn = nb d'outliers gènes /outsp = nb d'ouliers sp
##sp= 0 -> HGT sur tous les branches / sp = 1 -> HGT que sur les branches extérieures
SimOutliersLg <-function(nbgn, nbsp, outgn, outsp, outcell, sp = 1){
  tree=rtree(nbsp,rooted = TRUE, min=1,max=10)
  write.tree(tree, file = "arbre.tree")
  ListOutGnTree =list()
  "multiPhylo"->class(ListOutGnTree)
  for (i in 1:nbgn){
    ListOutGnTree[[i]]=tree
  }
  if(outsp!=0){ 
    ## outspe = seulement les branches externes
    if(sp==1){
      samp = sample(tree$tip.label,outsp)
      print(paste("outsp",samp,sep=" "))
      species=samp
      for (i in 1:length(samp)){
        j<-which(tree$tip.label==samp[i])
        l<-which(tree$edge[,2]==j)
        ListOutGnTree = BrLengthSp(ListOutGnTree, l)
      }
    }
    ## outspe = toutes les branches sont considérées et pas seulement les branches externes
    else if(sp==0){
      samp = sample(1:nrow(tree$edge),outsp)
      print(paste("outsp",samp,sep=" "))
      species=samp
      for (s in 1:outsp){
        br<-tree$edge[samp[s],2]
        print(tree$tip.label[Descendants(tree, br, "tips")[[1]]])
        ListOutGnTree =  BrLengthSp(ListOutGnTree, samp[s])
      }
    }
  }
  else {species = NULL}
  if (outgn !=0){
    samp= sample(1:length(ListOutGnTree),outgn)
    print(paste("outgn",samp,sep=" "))
    genes = samp
    for (j in 1:outgn){
      ListOutGnTree[[samp[[j]]]] =  BrLGn(ListOutGnTree[[samp[[j]]]], b=nrow(ListOutGnTree[[samp[[j]]]]$edge))
    }
  }
  else {genes = NULL}
  if (outcell !=0){## Les outcell ne peuvent pas être sur les espèces ou sur les gènes déjà complete outlier. 
    #Si on veut que ce soit possible d'avoir des outcell sur des gènes ou sp outliers complets, on fixe outcell = 0 et on créé des outcells avec la fonction
    #BrLengthOutCell()
    if (!is.null(genes)){genespossible = c(1:length(ListOutGnTree))[which(c(1:length(ListOutGnTree))!=genes)]}
    else{genespossible=c(1:length(ListOutGnTree))}
    samp= sample(genespossible,outcell)
    for (T in 1:outcell){
      if (!is.null(species)){speciespossible = setdiff(ListOutGnTree[[samp[T]]]$tip.label, species)}
      else{speciespossible=ListOutGnTree[[samp[T]]]$tip.label}
      lab = sample(speciespossible,1)
      print(paste("outcell =", samp[T],"/", lab[[1]],sep=" "))
      j<-which(ListOutGnTree[[samp[T]]]$tip.label==lab[[1]])
      l<-which(ListOutGnTree[[samp[T]]]$edge[,2]==j)
      ratio = runif(1,2,10)
      ListOutGnTree[[samp[T]]] = BrLength(ListOutGnTree[[samp[T]]], l ,ratio)
    }
  }
  ListTreesOut=list()
  "multiPhylo"->class(ListTreesOut)
  compteur=0
  for  (i in 1: length(ListOutGnTree)){
    tree<-di2multi(ListOutGnTree[[i]], tol=1e-3)
    n=length(tree$tip.label)
    write.tree(tree, file = "arbreHGT.tree")
    system(paste("perl WriteRose.pl -a arbreHGT.tree -p roseParam -l ", n, sep=""))
    system("rose param.output")
    system("phyml -i RoseTree.phy -m JC69 -n 1 -o tlr --quiet > sortie.out")
    ListTreesOut[[i]]= read.tree(file="RoseTree.phy_phyml_tree")
    compteur=compteur+1
    print(compteur)
  }
  return(ListTreesOut)
}

#--------------------------HGT----------------------------------------------------

##Fonction qui effectue n transferts horizontaux aléatoires dans un arbre de gène --> outgn
HGT1gn <- function(Tree, n=30){
  Tree2=Tree
  for (i in 1:n){
    Tree2 = HGT2(Tree2)
  }
  return(Tree2)
}

##--> créé k OutCell dans une liste d'arbres
HGToutCell <- function(ListTrees, k=1){
  ListTrees2 = ListTrees
  samp= sample(1:length(ListTrees),k)
  for (i in 1:k){
    lab = sample(ListTrees[[samp[i]]]$tip.label,1)
    print(samp[i])
    print(lab[[1]])
    ListTrees2[[samp[i]]] = HGT2(ListTrees[[samp[i]]],lab[[1]])
  }
  return(ListTrees2)
}

##Fonction qui créés k outliers gènes dans une liste d'arbres --> outgn
##n est le nombre de HGT par arbre (ne doit pas dépasser le nombre de branche de l'arbre)
HGToutgn <- function(ListTrees, n=30, k=1){
  ListTrees2 = ListTrees
  samp = sample(1:length(ListTrees),k)
  print(samp)
  for (j in 1:length(samp)){
    ListTrees2[[samp[j]]] =  HGT1gn(ListTrees2[[samp[j]]], n)
  }
  return(ListTrees2)
}

##Fonction qui change dans tous les arbres la même espèce (ou le groupe d'espèces selon la branche choisie) mais pas au même endroit --> outsp
HGToutsp <- function(ListTrees, species){ #ne pas laisser species = NULL sinon on a pas toujours la même branche qui bouge
  ListTrees2=ListTrees
  for (T in 1:length(ListTrees)){
    treeHGT = HGT2(ListTrees[[T]], species)
    ListTrees2[[T]]=treeHGT
  }
  return(ListTrees2)
}

## Fonction qui effectue un trasfert horizontal d'un endroit de l'arbre à un autre a partir d'une espèce ou d'un groupe d'espèces qu'on lui donne en entrée.
## Dans le cas où aucune espèce n'est précisée, le moment de la coupure est aléatoire.
## La branche coupée peut-être rebranché sur n'importe quelle branche entre le temps 0 et le temps de la coupure.
HGT2 <-  function(Tree, species=NULL){
  nbSpTot = Ntip(Tree)
  distnodes<-dist.nodes(Tree)[,nbSpTot+1] ##on ne le calcule qu'une seule fois et on ne garde que la colonne intéressante (distance à la racine)
  matricePos = cbind(distnodes[Tree$edge[,1]], distnodes[Tree$edge[,2]]) 
  rownames(matricePos)<-NULL
  if (!is.null(species)){
    if (length(species)==1) branche<-which(Tree$edge[,2]==which(Tree$tip.label==species))
    else branche<-which(Tree$edge[,2]==getMRCA(Tree, species))
    out = as.integer(branche)
    noeudOut=Tree$edge[out,2] ## noeud donneur
    ##moment de la coupure
    N <- runif(1, min = matricePos[branche,1], max= matricePos[branche,2]) #Choix d'un temps N aléatoire sur la branche choisie
    ##Branches coupées à N:
    ListNumBranche<-which(((matricePos[,1]<=N)&(matricePos[,2]>=N)))
    ListNumBranche2<-which(((matricePos[,1]<=N)))
  }
  #Si aucune espèce ou liste d'espèces n'est proposée, la coupure se fait aléatoirement dans le temps
  else{
    ##moment de la coupure
    N <- runif(1, min = min(distnodes), max= max(distnodes))
    ##Branches coupées à N:
    ListNumBranche<-which(((matricePos[,1]<=N)&(matricePos[,2]>=N)))
    ListNumBranche2<-which(((matricePos[,1]<=N)))
    if (length(ListNumBranche)==1) {
      branche = ListNumBranche[[1]]
      out = as.integer(branche)
      noeudOut=Tree$edge[out,2] ## noeud donneur
    }
    else {
      branche = sample(ListNumBranche,1)
      out = as.integer(branche)
      noeudOut=Tree$edge[out,2] ## noeud donneur
    }
  }
  if (length(ListNumBranche2)<=1) {
    return(Tree)
  }
  else if (length(ListNumBranche2)==2) {
    ins = as.integer(ListNumBranche2[ListNumBranche2!=as.integer(out)])
    noeudIns=Tree$edge[ins,2] ## noeud receveur
  }
  else if (length(ListNumBranche2)>2){
    ins = as.integer(sample(ListNumBranche2[ListNumBranche2!=as.integer(out)],1))
    noeudIns=Tree$edge[ins,2] ## noeud receveur
  }   
  ##création du sous-arbre qui va bouger
  b=rtree(2)
  tiplab<-c("branch", "out")
  b$tip.label=tiplab
  Min = vector()
  if (N > matricePos[ins,1]){Min[length(Min)+1]=N} 
  if (matricePos[ins,2] > matricePos[ins,1]){Min[length(Min)+1]=matricePos[ins,2]} 
  if (matricePos[out,2] > matricePos[ins,1]){Min[length(Min)+1]=matricePos[out,2]} 
  N2=runif(1, min = matricePos[ins,1], max= min(Min))
  ##Si la branche coupée est une banche interne
  if(noeudOut>nbSpTot){
    b$edge.length=c(dist.nodes(Tree)[noeudOut,nbSpTot+1]-N2,0) #Longueur de la branche choisie après avoir été coupée
    treeex <- extract.clade(Tree,noeudOut)
    bindedBranch = bind.tree(b, treeex, where = 1)
    Names = list()
    for (i in 1: length(Tree$tip.label)){
      if (is.element(Tree$tip.label[i], treeex$tip.label)){
        Names[i] <- "ancienNoeud"
      }
      else {
        Names[i] <- Tree$tip.label[i]
      }
    }
  }
  ##Si la branche coupée est une banche contenant une feuille à son extrémité (branche terminale)
  if(noeudOut<=nbSpTot){
    b$edge.length=c(0,0)
    c=rtree(2)
    tiplab2<-c(Tree$tip.label[noeudOut], "out")
    c$tip.label=tiplab2
    c$edge.length=c(dist.nodes(Tree)[noeudOut,nbSpTot+1]-N2,0) #Longueur de la branche choisie après avoir été coupée
    bindedBranch = bind.tree(b, c, where = 1)
    Names = list()
    for (i in 1: length(Tree$tip.label)){
      if (Tree$tip.label[[i]]==Tree$tip.label[[noeudOut]]){
        Names[i] <- "ancienNoeud"
      }
      else {
        Names[i] <- Tree$tip.label[i]
      }
    }
  }
  #On ne modifie pas l'arbre de départ
  tree3<-Tree
  tree3$tip.label<-Names
  treebinded <- bind.tree(tree3, bindedBranch, where=noeudIns, position = dist.nodes(Tree)[noeudIns,nbSpTot+1]-N2)
  treebinded <- drop.tip(treebinded,"out")
  treebinded <- drop.tip(treebinded,"ancienNoeud")
  return(treebinded)
}

#--------------------------Longueur de Branche---------------------------------------------------

##Fonction qui permet de changer la longueur d'une branche donnée dans un arbre donné en multipliant sa longeur par un ratio
BrLength <- function(Tree, branche = sample(1:nrow(Tree$edge),1), ratio){
  Tree2 = Tree
  Tree2$edge.length[branche]=Tree$edge.length[branche]*ratio
  return(Tree2)
}

##Fonction qui permet d'insérer des outliers cells dans une liste d'arbres
BrLengthOutCell <- function(ListTrees, k=1, ratio){
  ListTrees2=ListTrees
  samp = sample(1:length(ListTrees),k)
  for (T in 1:k){
    lab = sample(ListTrees[[samp[T]]]$tip.label,1)
    print(samp[T])
    print(lab)
    j<-which(ListTrees[[samp[T]]]$tip.label==lab)
    l<-which(ListTrees[[samp[T]]]$edge[,2]==j)
    ListTrees2[[samp[T]]] = BrLength(ListTrees[[samp[T]]], l ,ratio)
  }
  return(ListTrees2)
}

##Fonction qui permet de changer aléatoirement la longueur d'une branche donnée dans tous les arbres d'une liste -> outsp
BrLengthSp <- function(ListTrees, branche){
  ListTrees2=ListTrees
  for (t in 1:length(ListTrees)){
    Tree2 = ListTrees[[t]]
    ratio = runif(1, 1.5, 3)
    Tree2 = BrLength(Tree2, branche, ratio)
    ListTrees2[[t]]=Tree2
  }
  return(ListTrees2)
}

##Fonction qui permet de changer aléatoirement la longueur d'un certain nombre (b) de branches d'un arbre -> outgn
BrLGn <- function(Tree, b){
  Tree2=Tree
  branche = sample(1:nrow(Tree$edge),b)
  for (i in 1:b){
    ratio = runif(1, 1.5, 3)
    Tree2 = BrLength(Tree2, branche[[i]], ratio)
  }
  return(Tree2)
}

##Changer la longueur d'un certain nombre de branche de plusieurs arbres dans une liste
##b = nombre de branches modifiées dans un arbre, k = nombre d'arbres modifiés
BrLengthGn <- function(ListTrees, b, k){
  ListTrees2 = ListTrees
  samp= sample(1:length(ListTrees),k)
  print(samp)
  for (j in 1:k){
    ListTrees2[[samp[[j]]]] =  BrLGn(ListTrees[[samp[[j]]]], b)
  }
  return(ListTrees2)
}


