
######-----------------------------------Simulation--------------------------------------------------------------------------

###Fonction qui génère une liste d'arbres de nbgn gènes avec nbsp espèces contenant des outliers gènes (outgn) et espèces (outsp) générés par HGT
##nbsp = nombre d'espèces dans l'arbre / nbgn = nombre d'arbres / outgn = nb d'outlier gènes /outsp = nb d'oulier sp  /outcell = couple outlier gn/sp
##sp= 0 -> HGT sur tous les branches / sp = 1 -> HGT que sur les branches extérieures
SimOutliersHGT <-function(tree, nbgn, outgn, outsp, outcell, sp = 0){
  genes=NULL
  nbsp = length(tree$tip.label)
  #tree<-rtree(nbsp,rooted = TRUE,min=1,max=10)
  write.tree(tree, file = "arbre.tree")
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
        #system(paste("echo", paste("outsp",samp,sep=" "), ">>log_",i_plan, sep=" "))
        species = samp
        ListOutGnTree = HGToutsp(ListOutGnTree,samp)
      }
    }
  }
  else {species = NULL}
  if (outgn !=0){
    samp = sample(1:length(ListOutGnTree),outgn)
    print(paste("outgn",samp,sep=" "))
    #system(paste("echo", paste("outgn",samp,sep=" "), ">>log_",i_plan,sep=" "))
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
      #system(paste("echo", paste("outcell =", samp[i],"/", lab[[1]],sep=" "), ">>log_" ,i_plan, sep=" "))
      ListOutGnTree[[samp[i]]] = HGT2(ListOutGnTree[[samp[i]]],lab[[1]])
    }
  }
  ListTreesOut=list()
  "multiPhylo"->class(ListTreesOut)
  #compteur= 0
  for  (i in 1: length(ListOutGnTree)){
    tree<-di2multi(ListOutGnTree[[i]], tol=1e-3)
    n=length(tree$tip.label)
    write.tree(tree, file = "arbreHGT.tree")
    system(paste("perl WriteRose.pl -a arbreHGT.tree -p roseParam -l ", n, sep=""))
    system("rose param.output")
    system("phyml -i RoseTree.phy -m JC69 -n 1 -o tlr --quiet > sortie.out")
    #system("phyml -i RoseTree.phy -m JC69 -n 1 -o lr -u arbreHGT.tree --quiet > sortie.out")
    ListTreesOut[[i]]= read.tree(file="RoseTree.phy_phyml_tree")
    #compteur=compteur+1
    #print(compteur)
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
SimOutliersLg <-function(tree, nbgn, outgn, outsp, outcell, sp = 1){
  nbsp = length(tree$tip.label)
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
