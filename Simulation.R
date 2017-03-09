##load required packages
require(ape)
require(phangorn)

###Fonction qui génère une liste d'arbres de nbgn gènes avec nbsp espèces contenant des outliers gènes (outgn) et espèces (outsp) générés par HGT
##nbsp = nombre d'espèces dans l'arbre / nbgn = nombre d'arbres / outgn = nb d'outlier gènes /outsp = nb d'oulier sp 
SimOutliersHGT <-function(nbsp=10, nbgn=10, outgn=1, outsp=1, sp = NULL){
    tree<-rtree(nbsp,rooted = FALSE,min=2,max=5)
    write.tree(tree, file = "arbre.tree")
    ListOutGnTree =list()
    "multiPhylo"->class(ListOutGnTree)
    for (i in 1:nbgn){ # n fois le même arbre qui va evoluer différement
        ListOutGnTree[[i]]=tree
    }
    if(outsp!=0){
      if(!is.null(sp)){
        ## On ne fait des hgt que sur les branches externes (pour pouvoir contrôler le nombre d'espèces outliers)
        samp = sample(tree$tip.label,outsp) 
        print(samp)
        for (o in 1:outsp){
          ListOutGnTree = HGToutsp(ListOutGnTree, samp[o])
        }
      }
      else{
        s=1
        #dans le cas où outspe > 1, on peut prendre plusieurs fois la même branche.
        #C'est également possible que deux branches différentes choisies comme outliers aient les mêmes espèces déscendantes
        bra = c(1:nrow(tree$edge))
        while (s <= outsp){
            br<-sample(bra,1) ##une branche au hasard 
            br<-tree$edge[br,2]
            samp<-tree$tip.label[Descendants(tree, br, "tips")[[1]]]
            print(samp)
            ListOutGnTree = HGToutsp(ListOutGnTree,samp)
            s=s+1
        }
      }
    }
    if (outgn !=0){
        ListOutGnTree = HGToutgn(ListOutGnTree,n=nrow(tree$edge),k=outgn) #autant d'hgt que de branches par arbre
    }
  
    ListTreesOut=list()
    "multiPhylo"->class(ListTreesOut)
    for  (i in 1: length(ListOutGnTree)){
      
      write.tree(ListOutGnTree[[i]], file = "arbreHGT.tree")
      
      system(paste("perl WriteRose.pl -a arbreHGT.tree -p roseParam -l ", nbsp, sep=""))
      system("rose param.output")
      system("phyml -i RoseTree.phy -n 1 -o lr -u arbre.tree --quiet")
      ListTreesOut[[i]]= read.tree(file="RoseTree.phy_phyml_tree")
      
      #system("/home/aurore/Téléchargements/Seq-Gen.v1.3.3/source/seq-gen -mHKY85 -n1 -l100 < arbreHGT.tree > seqtrees.dat")
      #system("phyml -i seqtrees.dat -n 1 -o lr -u arbre.tree --quiet")
      #ListTreesOut[[i]]= read.tree(file="seqtrees.dat_phyml_tree")
      
    }
    #return(ListOutGnTree)
    return(ListTreesOut)
}

###Fonction qui génère une liste d'arbres de nbgn gènes avec nbsp espèces contenant des outliers gènes (outgn) et espèces (outsp) générés
###en modifiant les longueurs de branches
#nbsp = nombre d'espèces dans l'arbre / nbgn = nombre d'arbres / outgn = nb d'outliers gènes /outsp = nb d'ouliers sp 
SimOutliersLg <-function(nbsp, nbgn, outsp, outgn, sp = NULL){
  tree=rtree(nbsp,rooted = TRUE, min=2,max=5)
  write.tree(tree, file = "arbre.tree")
  ListOutGnTree =list()
  "multiPhylo"->class(ListOutGnTree)
  for (i in 1:nbgn){
    ListOutGnTree[[i]]=tree
  }
  if(outsp!=0){ 
    ## outspe = seulement les branches externes
    if(!is.null(sp)){
      samp = sample(tree$tip.label,outsp)
      print(samp)
      for (i in 1:length(samp)){
        j<-which(tree$tip.label==samp[i])
        l<-which(tree$edge[,2]==j)
        ListOutGnTree = BrLengthSp(ListOutGnTree, l)
      }
    }
    ## outspe = toutes les branches sont considérées et pas seulement les branches externes
    else{
      samp = sample(1:nrow(tree$edge),outsp)
      for (s in 1:outsp){
        br<-tree$edge[samp[s],2]
        print(tree$tip.label[Descendants(tree, br, "tips")[[1]]])
        ListOutGnTree =  BrLengthSp(ListOutGnTree, samp[s])
      }
    }
  }
  if (outgn !=0){
    ListOutGnTree = BrLengthGn(ListOutGnTree, b=nrow(tree$edge), k=outgn)
  }
  
  ListTreesOut=list()
  "multiPhylo"->class(ListTreesOut)
  for  (i in 1: length(ListOutGnTree)){
    
    write.tree(ListOutGnTree[[i]], file = "arbreHGT.tree")
    
    system(paste("perl WriteRose.pl -a arbreHGT.tree -p roseParam -l ", nbsp, sep=""))
    system("rose param.output")
    system("phyml -i RoseTree.phy -n 1 -o lr -u arbre.tree --quiet")
    ListTreesOut[[i]]= read.tree(file="RoseTree.phy_phyml_tree")
    
    #system("/home/aurore/Téléchargements/Seq-Gen.v1.3.3/source/seq-gen -mHKY85 -n1 -l100 < arbreHGT.tree > seqtrees.dat")
    #system("phyml -i seqtrees.dat -n 1 -o lr -u arbre.tree --quiet")
    #ListTreesOut[[i]]= read.tree(file="seqtrees.dat_phyml_tree")
    
  }
  #return(ListOutGnTree)
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

## NOUVELLE TENTATIVE. Je ne modifie pas ta fonction pour ne pas trop mettre le bazarre.
## je propose que la fonction fasse les transferts en fonction d'une espèce ou d'un groupe d'espèces
## plutôt qu'a partir d'un numéro de branche 
HGT2 <-  function(Tree, species=NULL){
  nbSpTot = Ntip(Tree)
  distnodes<-dist.nodes(Tree)[,nbSpTot+1] ##on ne le calcule qu'une seule fois et on ne garde que la colonne intéressante (distance à la racine)
  matricePos = cbind(distnodes[Tree$edge[,1]], distnodes[Tree$edge[,2]]) 
  rownames(matricePos)<-NULL
  if (!is.null(species)){
    if (length(species)==1) branche<-which(Tree$edge[,2]==which(Tree$tip.label==species))
    else branche<-which(Tree$edge[,2]==getMRCA(Tree, species))
    ##moment de la coupure
    N <- runif(1, min = matricePos[branche,1], max= matricePos[branche,2]) #Choix d'un temps N aléatoire sur la branche choisie
    ##Branches coupées à N:
    ListNumBranche<-which(((matricePos[,1]<=N)&(matricePos[,2]>=N)))
  }
  #Si aucune espèce ou liste d'espèces n'est proposée, la coupure se fait aléatoirement dans le temps
  else{
    ##moment de la coupure
    N <- runif(1, min = min(distnodes), max= max(distnodes))
    ##Branches coupées à N:
      ListNumBranche<-which(((matricePos[,1]<=N)&(matricePos[,2]>=N)))
      branche = as.integer(sample(ListNumBranche,1))
      
  }
  if (length(ListNumBranche)==1) {
      return(Tree)
  }
  else {  
    if (length(ListNumBranche)==2) {
      out = as.integer(branche)
      noeudOut=Tree$edge[out,2] ## noeud donneur
      ins = as.integer(ListNumBranche[ListNumBranche!=as.integer(out)])
      noeudIns=Tree$edge[ins,2] ## noeud receveur
    }
    else if (length(ListNumBranche)>2){
        out = as.integer(branche)
        noeudOut=Tree$edge[out,2] ## noeud donneur
        ins = as.integer(sample(ListNumBranche[ListNumBranche!=as.integer(out)],1))
        noeudIns=Tree$edge[ins,2] ## noeud receveur
    }   
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
        ##Si la branche coupée est une banche contenant une feuille à son extrémité (branche terminale)
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
        tree3<-Tree
        tree3$tip.label<-Names
        treebinded <- bind.tree(tree3, bindedBranch, where=noeudIns, position = dist.nodes(Tree)[noeudIns,nbSpTot+1]-N)
        treebinded <- drop.tip(treebinded,"out")
        treebinded <- drop.tip(treebinded,"ancienNoeud")
        return(treebinded)
    }
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
    ratiomin = runif(1, 0.1, 0.4)
    ratiomax = runif(1, 1.6, 1.9)
    ratio = sample(c(ratiomin, ratiomax),1)
    Tree2 = BrLength(Tree2, branche, ratio)
    ListTrees2[[t]]=Tree2
  }
  return(ListTrees2)
}

##Fonction qui permet de changer aléatoirement la longueur d'un certain nombre (b) de branches d'un arbre -> outgn
BrLGn <- function(Tree, b){
  Tree2=Tree
  for (i in 1:b){
    branche = sample(1:nrow(Tree$edge),1)
    ratiomin = runif(1, 0.4, 0.7)
    ratiomax = runif(1, 1.3, 1.6)
    ratio = sample(c(ratiomin, ratiomax),1)
    Tree2 = BrLength(Tree2, branche, ratio)
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
