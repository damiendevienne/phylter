#--------------------------HGT----------------------------------------------------

##Fonction qui éfectue n transfert horizontaux aléatoire dans un arbre de gène --> outgn
HGToutgn <- function(Tree, n=10){
  Tree2 = Tree
  for (i in 1:n){
    Tree2 = HGT(Tree2,branche = sample(1:nrow(Tree2$edge),1))
  }
  return(Tree2)
}

##Fonction qui change dans tous les arbres la même espèce (ou le groupe d'espèce selon la branche choisie) mais pas au même endroit --> outsp
HGToutsp <- function(ListTrees, branche){
  ListTree2=list()
  for (t in 1:length(ListTrees)){
    treeHGT = HGT(ListTrees[[T]], branche)
    ListTree2[[T]]=treeHGT
  }
  return(ListTrees)
}

##Cette fonction permet d'effectuer un transfert dans un arbre en choisissant sur quelle branche se fait la coupure au départ. 
##(numéro de la branche = numéro de la ligne dans tree$edge).
##Le branchement du sous-arbre coupé se fera sur une autre branche (mais au même temps relatif)
HGT <-  function(Tree, branche = sample(1:nrow(Tree$edge),1)){
    nbSpTot = length(Tree$tip.label)
    matricePos = matrix(nrow=nrow(Tree$edge), ncol=ncol(Tree$edge))
    for (i in 1:nrow(Tree$edge)){
      for (j in 1:ncol(Tree$edge)){
        matricePos[i,j]=dist.nodes(Tree)[Tree$edge[i,j],nbSpTot+1]
      }
    }
    #Il faut qu'il y ai au moins deux point de coupure à un temps t pour réaliser un déplacement (un acctepeur et un donneur):
    #Cela peut être une raison pour que la fonction bugge
    ListNumBranche <- list()
    while (length(ListNumBranche) <= 1){
      N <- runif(1, min = matricePos[branche,1], max= matricePos[branche,2]) #Choix d'un temps N aléatoire sur la branche choisie
      for (i in 1:nrow(matricePos)){
        if((matricePos[i,1] < N) & (matricePos[i,2] > N)){
          ListNumBranche[length(ListNumBranche)+1] <- i #Est ce qu'une autre branche dans l'arbre comprend ce temps N
        }
      }
    }
      out = as.numeric(sample(ListNumBranche,1)) 
      noeudOut=Tree$edge[out,2] ## noeud donneur
      ins = as.numeric(sample(ListNumBranche[ListNumBranche!=as.numeric(out)],1))
      noeudIns=Tree$edge[ins,2] ## noeud receveur
      ##création du sous-arbre qui va bouger
      b=rtree(2)
      tiplab<-c("branch", "out")
      b$tip.label=tiplab
      ##Si la branche coupée est une banche interne
      if(noeudOut>nbSpTot){
        b$edge.length=c(dist.nodes(Tree)[noeudOut,nbSpTot+1]-N,0) #Longueur de la branche choisie après avoir été coupée
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
      ##Si la branche coupée est une banche contenant une feuille à son extrémité
      if(noeudOut<=nbSpTot){
        b$edge.length=c(0,0)
        c=rtree(2)
        tiplab2<-c(Tree$tip.label[noeudOut], "out")
        c$tip.label=tiplab2
        c$edge.length=c(dist.nodes(Tree)[noeudOut,nbSpTot+1]-N,0) #Longueur de la branche choisie après avoir été coupée
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
      tree3=Tree
      tree3$tip.label=Names
      treebinded = bind.tree(tree3, bindedBranch, where=noeudIns, position = dist.nodes(Tree)[noeudIns,nbSpTot+1]-N)
      treebinded = drop.tip(treebinded,"out")
      treebinded = drop.tip(treebinded,"ancienNoeud")
      return(treebinded)
}

#--------------------------Longueur de Branche---------------------------------------------------

##Fonction qui permet de changer la longueur d'une branche donnée dans un arbre donné
BrLength <- function(Tree, branche = sample(1:nrow(Tree$edge),1), ratio = 2){
  Tree2 = Tree
  Tree2$edge.length[branche]=Tree2$edge.length[branche]*ratio
  return(Tree2)
}

##Fonction qui permet de changer aléatoirement la longueur d'une branche donnée dans tous les arbres d'une liste -> outsp
BrLengthSp <- function(ListTrees, branche, ratioMin = 0.1, ratioMax = 10){
  ListTrees2=list()
  for (t in 1:length(ListTrees)){
    ratio = runif(1,ratioMin,ratioMax)
    Tree2 = BrLength(ListTrees[[t]], branche, ratio)
    ListTrees2[[t]]=Tree2
  }
  return(ListTrees2)
}

##Fonction qui permet de changer aléatoirement la longueur d'un certain nombre de branches d'un arbre -> outgn
BrLengthGn <- function(Tree, n=10, ratioMin = 0.1, ratioMax = 10){
  Tree2=Tree
  for (i in 1:n){
    branche = sample(1:nrow(Tree2$edge),1)
    ratio = runif(1,ratioMin,ratioMax)
    Tree2 = BrLength(Tree2, branche, ratio)
  }
  return(Tree2)
}


