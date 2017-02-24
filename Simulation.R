###Fonction qui génère une liste d'arbre de nbgn gènes avec nbsp espèces contenant des outliers gènes (Outgn) et espèces (Outsp) générés
###en générant des HGT
SimOutliersHGT <-function(nbsp = 30, nbgn = 30, OutSp = 1, Outgn= 1){
  tree=rtree(nbsp,rooted = TRUE)
  write.tree(tree, file = "arbre.tree")
  ListOutGnTree =list()
  for (i in 1:nbgn){
    ListOutGnTree[[i]]=tree
  }
  if(OutSp!=0){
    ## outspe = seulement les feuilles
    samp = sample(1:nbsp,OutSp)
    for (s in 1:length(samp)){
      for (j in 1:nrow(tree$edge)){ 
        if (tree$edge[j,2]==samp[[s]]){
          ListOutGnTree = HGToutsp(ListOutGnTree, j)
        }
      }
    }
  }
  if (Outgn !=0){
    samp = sample(1:nbgn,Outgn)
    for (i in 1:length(samp)){
      ListOutGnTree[[samp[i]]] = HGT1gn(ListOutGnTree[[samp[i]]], n=30)
    }
  }
  ListTreesOut=list()
  for  (i in 1: length(ListOutGnTree)){
    write.tree(ListOutGnTree[[i]], file = "arbreHGT.tree")
    system("/home/aurore/Téléchargements/Seq-Gen.v1.3.3/source/seq-gen -mGTR -n1 -s0.5 < arbreHGT.tree > seqtrees.dat")
    system("phyml -i seqtrees.dat -n 1 -o lr -u arbre.tree --quiet")
    ListTreesOut[[i]]= read.tree(file="seqtrees.dat_phyml_tree")
  }
  return(ListTreesOut)
}

###Fonction qui génère une liste d'arbre de nbgn gènes avec nbsp espèces contenant des outliers gènes (Outgn) et espèces (Outsp) générés
###en modifiant les longueurs de branches
SimOutliersLg <-function(nbsp = 30, nbgn = 30, OutSp = 1, Outgn= 1){
  tree=rtree(nbsp,rooted = TRUE)
  write.tree(tree, file = "arbre.tree")
  ListOutGnTree =list()
  for (i in 1:nbgn){
    ListOutGnTree[[i]]=tree
  }
  if(OutSp!=0){
    ## outspe = seulement les feuilles
    samp = sample(1:nbsp,OutSp)
    for (s in 1:length(samp)){
      for (j in 1:nrow(tree$edge)){ 
        if (tree$edge[j,2]==samp[[s]]){
          ListOutGnTree = BrLengthSp(ListOutGnTree, j, ratiomin = 0.1, ratiomax=10)
        }
      }
    }
  }
  if (Outgn !=0){
    samp = sample(1:nbgn,Outgn)
    ListOutGnTree = BrLengthGn(ListOutGnTree,n=30,Listgn=samp, ratio=3)
  }
  ListTreesOut=list()
  for  (i in 1: length(ListOutGnTree)){
    write.tree(ListOutGnTree[[i]], file = "arbreHGT.tree")
    system("/home/aurore/Téléchargements/Seq-Gen.v1.3.3/source/seq-gen -mGTR -n1 -s0.5 < arbreHGT.tree > seqtrees.dat")
    system("phyml -i seqtrees.dat -n 1 -o lr -u arbre.tree --quiet")
    ListTreesOut[[i]]= read.tree(file="seqtrees.dat_phyml_tree")
  }
  return(ListTreesOut)
}


#--------------------------HGT----------------------------------------------------

##Fonction qui effectue n transfert horizontaux aléatoire dans un arbre de gène --> outgn
HGT1gn <- function(Tree, n=15){
  Tree2 = Tree
  for (i in 1:n){
    Tree2 = HGT(Tree,branche = sample(1:nrow(Tree$edge),1))
  }
  return(Tree2)
}

##--> OutCell
HGToutCell <- function(ListTrees, n=1){
  ListTrees2 = ListTrees
  samp= sample(1:length(ListTrees),n)
  for (i in 1:n){
    ListTrees2[[samp[[i]]]] = HGT(ListTrees[[samp[[i]]]])
  }
  return(ListTrees2)
}

##Fonction qui créés des outliers gènes dans une liste d'arbre --> outgn
HGToutgn <- function(ListTrees, n=15, k=1, Listgn = NULL){
  ListTrees2 = ListTrees
  if (is.null(Listgn)){
    samp= sample(1:length(ListTrees),k)
    for (j in 1:length(samp)){
      ListTrees2[[samp[[j]]]] =  HGT1gn(ListTrees[[samp[[j]]]], n)
    }
  }
  else{
    for (i in 1:length(Listgn)){
      ListTrees2[[Listgn[[i]]]] =  HGT1gn(ListTrees[[Listgn[[i]]]], n)
    }
  }
  return(ListTrees2)
}

##Fonction qui change dans tous les arbres la même espèce (ou le groupe d'espèce selon la branche choisie) mais pas au même endroit --> outsp
HGToutsp <- function(ListTrees, branche){
  ListTrees2=ListTrees
  for (T in 1:length(ListTrees)){
    treeHGT = HGT(ListTrees[[T]], branche)
    ListTrees2[[T]]=treeHGT
  }
  return(ListTrees2)
}

##Cette fonction permet d'effectuer un transfert dans un arbre en choisissant sur quelle branche se fait la coupure au départ. 
##(numéro de la branche = numéro de la ligne dans tree$edge).
##Le branchement du sous-arbre coupé se fera sur une autre branche (mais au même temps relatif)
HGT <-  function(Tree, branche = sample(1:nrow(Tree$edge),1)){
  branche = sample(1:nrow(Tree$edge),1)
    nbSpTot = length(Tree$tip.label)
    matricePos = matrix(nrow=nrow(Tree$edge), ncol=ncol(Tree$edge))
    for (i in 1:nrow(Tree$edge)){
      for (j in 1:ncol(Tree$edge)){
        matricePos[i,j]=dist.nodes(Tree)[Tree$edge[i,j],nbSpTot+1]
      }
    }
    #Il faut qu'il y ai au moins deux point de coupure à un temps t pour réaliser un déplacement (un acctepeur et un donneur):
    Boo = FALSE
    while (Boo == FALSE){
      ListNumBranche <- list()
      N <- runif(1, min = matricePos[branche,1], max= matricePos[branche,2]) #Choix d'un temps N aléatoire sur la branche choisie
      for (i in 1:nrow(matricePos)){
        if((matricePos[i,1] < N) & (matricePos[i,2] > N) & (!(i%in%ListNumBranche))){
          ListNumBranche[length(ListNumBranche)+1] <- i #Est ce qu'une autre branche dans l'arbre comprend ce temps N
        }
      }
      if (length(ListNumBranche) > 1){
        Boo = TRUE
      }
    }  
      out = as.integer(sample(ListNumBranche,1)) 
      noeudOut=Tree$edge[out,2] ## noeud donneur
      ins = as.integer(sample(ListNumBranche[ListNumBranche!=as.integer(out)],1))
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
BrLength <- function(Tree, branche = sample(1:nrow(Tree$edge),1), ratio=5){
  Tree2 = Tree
  Tree2$edge.length[branche]=Tree2$edge.length[branche]*ratio
  return(Tree2)
}

##Fonction qui permet d'insérer des outliers cells dans une liste d'arbres
BrLengthOutCell <- function(ListTrees, n=1, ratio=5){
  ListTrees2=ListTrees
  samp= sample(1:length(ListTrees),n)
  for (T in 1:n){
      ListTrees2[[samp[[T]]]] = BrLength(ListTrees[[samp[[T]]]],ratio)
  }
  return(ListTrees2)
}

##Fonction qui permet de changer aléatoirement la longueur d'une branche donnée dans tous les arbres d'une liste -> outsp
BrLengthSp <- function(ListTrees, branche, ratiomin = 0.1, ratiomax=10){
  ListTrees2=ListTrees
  for (t in 1:length(ListTrees)){
    ratio=runif(1,ratiomin,ratiomax)
    Tree2 = BrLength(ListTrees[[t]], branche, ratio)
    ListTrees2[[t]]=Tree2
  }
  return(ListTrees2)
}

##Fonction qui permet de changer aléatoirement la longueur d'un certain nombre de branches d'un arbre -> outgn
BrLGn <- function(Tree, n=20, ratio=5){
  Tree2=Tree
  for (i in 1:n){
    branche = sample(1:nrow(Tree2$edge),1)
    Tree2 = BrLength(Tree2, branche, ratio)
  }
  return(Tree2)
}

##Changer la longueur d'un certain nombre de branche de plusieurs arbres dans une liste
BrLengthGn <- function(ListTrees, n=20, k, Listgn=NULL, ratio=5){
  ListTrees2 = ListTrees
  if (is.null(Listgn)){
      samp= sample(1:length(ListTrees),k)
      for (j in 1:length(samp)){
        ListTrees2[[samp[[j]]]] =  BrLGn(ListTrees[[samp[[j]]]], n, ratio)
      }
  }
  else{
    for (i in 1:length(Listgn)){
      ListTrees2[[Listgn[[i]]]] =  BrLGn(ListTrees[[Listgn[[i]]]],n, ratio)
    }
  }
  return(ListTrees2)
}



