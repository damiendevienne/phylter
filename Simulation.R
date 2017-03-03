###Fonction qui génère une liste d'arbres de nbgn gènes avec nbsp espèces contenant des outliers gènes (outgn) et espèces (outsp) générés par HGT
#nbsp = nombre d'espèces dans l'arbre / nbgn = nombre d'arbres / outgn = nb d'outlier gènes /outsp = nb d'oulier sp 
#sp = "f" si on ne veut que les HGT (sp) ne se fasse que sur les branches externes (sinon ne rien mettre)
SimOutliersHGT <-function(nbsp=10, nbgn=10, outgn=1, outsp=1){
  tree=rtree(nbsp,rooted = TRUE,min=1,max=5)
  write.tree(tree, file = "arbre.tree")
  ListOutGnTree =list()
  "multiPhylo"->class(ListOutGnTree)
  for (i in 1:nbgn){ # n fois le même arbre qui va evoluer différement
    ListOutGnTree[[i]]=tree
  }
  if(outsp!=0){
    ## On ne fait des hgt que sur les branches externes (pour pouvoir contrôler le nombre d'espèces outliers)
    s=1
    while (s <= outsp){
      samp = sample(tree$tip.label,1)
      for (j in 1:length(tree$tip.label)){
        if (tree$tip.label[j]==samp){
          for (l in 1:nrow(tree$edge)){
            if (tree$edge[l,2]==j){
              ListOutGnTree = HGToutsp(ListOutGnTree, l)
              s=s+1
            }
          }
        }
      }
    }
  }
  if (outgn !=0){
    ListOutGnTree = HGToutgn(ListOutGnTree,n=nrow(tree$edge),k=outgn) #autant d'hgt que de branches par arbre
  }
  return(ListOutGnTree)
}

###Fonction qui génère une liste d'arbres de nbgn gènes avec nbsp espèces contenant des outliers gènes (outgn) et espèces (outsp) générés
###en modifiant les longueurs de branches
#nbsp = nombre d'espèces dans l'arbre / nbgn = nombre d'arbres / outgn = nb d'outliers gènes /outsp = nb d'ouliers sp 
SimOutliersLg <-function(nbsp, nbgn, outsp, outgn, sp="f"){
  tree=rtree(nbsp,rooted = TRUE, min=1,max=5)
  write.tree(tree, file = "arbre.tree")
  ListOutGnTree =list()
  "multiPhylo"->class(ListOutGnTree)
  for (i in 1:nbgn){
    ListOutGnTree[[i]]=tree
  }
  if(outsp!=0){ 
    ## outspe = seulement les branches externes
    if(sp == "f"){
      s=1
      while (s <= outsp){
        samp = sample(1:nbsp,1)
        for (j in 1:nrow(tree$edge)){ 
          if (tree$edge[j,2]==samp){
            ListOutGnTree = BrLengthSp(ListOutGnTree, j)
            s=s+1
          }
        }
      }
    }
    ## outspe = toutes les branches sont considérées et pas seulement les branches externes
    else{
      samp = sample(1:nrow(tree$edge),outsp)
      for (s in 1:outsp){
        ListOutGnTree =  BrLengthSp(ListOutGnTree, samp[[s]])
      }
    }
  }
  if (outgn !=0){
    ListOutGnTree = BrLengthGn(ListOutGnTree, b=nrow(tree$edge), k=outgn) #autant de changements de longueur que de branches par arbre
  }
  return(ListOutGnTree)
}

#--------------------------HGT----------------------------------------------------

##Fonction qui effectue n transferts horizontaux aléatoires dans un arbre de gène --> outgn
HGT1gn <- function(Tree, n=30){
  listBranche = sample(1:nrow(Tree$edge),n)
  for (i in 1:n){
    Tree2 = HGT(Tree,branche = listBranche[[i]])
  }
  return(Tree2)
}

##--> créé k OutCell dans une liste d'arbres
HGToutCell <- function(ListTrees, k=1){
  ListTrees2 = ListTrees
  samp= sample(1:length(ListTrees),k)
  for (i in 1:k){
    ListTrees2[[samp[[i]]]] = HGT(ListTrees[[samp[[i]]]])
  }
  return(ListTrees2)
}

##Fonction qui créés k outliers gènes dans une liste d'arbres --> outgn
##n est le nombre de HGT par arbre (ne doit pas dépasser le nombre de branche de l'arbre)
HGToutgn <- function(ListTrees, n=30, k=1){
  ListTrees2 = ListTrees
  samp= sample(1:length(ListTrees),k)
  for (j in 1:length(samp)){
     ListTrees2[[samp[[j]]]] =  HGT1gn(ListTrees2[[samp[[j]]]], n)
  }
  return(ListTrees2)
}

##Fonction qui change dans tous les arbres la même espèce (ou le groupe d'espèces selon la branche choisie) mais pas au même endroit --> outsp
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
    nbSpTot = length(Tree$tip.label)
    matricePos = matrix(nrow=nrow(Tree$edge), ncol=ncol(Tree$edge))
    for (i in 1:nrow(Tree$edge)){
      for (j in 1:ncol(Tree$edge)){
        matricePos[i,j]=dist.nodes(Tree)[Tree$edge[i,j],nbSpTot+1]
      }
    }
    #Il faut qu'il y ai au moins deux points de coupure à un temps t pour réaliser un déplacement (un accepteur et un donneur):
    Boo = FALSE
    while (Boo == FALSE){
      ListNumBranche <- list()
      N <- runif(1, min = matricePos[branche,1], max= matricePos[branche,2]) #Choix d'un temps N aléatoire sur la branche choisie
      for (i in 1:nrow(matricePos)){
        if((matricePos[i,1] < N) & (matricePos[i,2] > N) & (!(i%in%ListNumBranche))){
          ListNumBranche[length(ListNumBranche)+1] <- i #Est ce qu'une autre branche dans l'arbre comprend ce temps N?
        }
      }
      if (length(ListNumBranche) > 1){
        Boo = TRUE
      }
    }  
      out = as.integer(branche)
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

##Fonction qui permet de changer la longueur d'une branche donnée dans un arbre donné en multipliant sa longeur par un ratio
BrLength <- function(Tree, branche = sample(1:nrow(Tree$edge),1), ratio){
  Tree2 = Tree
  Tree2$edge.length[branche]=Tree$edge.length[branche]*ratio
  return(Tree2)
}

##Fonction qui permet d'insérer des outliers cells dans une liste d'arbres
BrLengthOutCell <- function(ListTrees, k=1, ratio){
  ListTrees2=ListTrees
  samp= sample(1:length(ListTrees),k)
  for (T in 1:k){
      ListTrees2[[samp[[T]]]] = BrLength(ListTrees[[samp[[T]]]], branche = sample(1:nrow(Tree$edge),1),ratio)
  }
  return(ListTrees2)
}

##Fonction qui permet de changer aléatoirement la longueur d'une branche donnée dans tous les arbres d'une liste -> outsp
BrLengthSp <- function(ListTrees, branche){
  ListTrees2=ListTrees
  for (t in 1:length(ListTrees)){
    ratiomin = runif(1, 0.3, 0.7)
    ratiomax = runif(1, 1.7, 2.3)
    ratio = sample(c(ratiomin, ratiomax),1)
    Tree2 = BrLength(ListTrees[[t]], branche, ratio)
    ListTrees2[[t]]=Tree2
  }
  return(ListTrees2)
}

##Fonction qui permet de changer aléatoirement la longueur d'un certain nombre (b) de branches d'un arbre -> outgn
BrLGn <- function(Tree, b){
  Tree2=Tree
  branche = sample(1:nrow(Tree$edge),b)
  for (i in 1:b){
    ratiomin = runif(1, 0.3, 0.7)
    ratiomax = runif(1, 1.7, 2.3)
    ratio = sample(c(ratiomin, ratiomax),1)
    Tree2 = BrLength(Tree2, branche[[i]], ratio)
  }
  return(Tree2)
}

##Changer la longueur d'un certain nombre de branche de plusieurs arbres dans une liste
##b = nombre de branches modifiées dans un arbre, k = nombre d'arbres modifiés
BrLengthGn <- function(ListTrees, b, k){
  ListTrees2 = ListTrees
  samp= sample(1:length(ListTrees),k)
  for (j in 1:k){
    ListTrees2[[samp[[j]]]] =  BrLGn(ListTrees[[samp[[j]]]], b)
  }
  return(ListTrees2)
}
