

#################1 transfert dans 1 arbre

##Choix de la branche

HGT <-  function(Tree, branche = sample(1:nrow(Tree$edge),1)){
    nbSpTot = length(Tree$tip.label)
    matricePos = matrix(nrow=nrow(Tree$edge), ncol=ncol(Tree$edge))
    for (i in 1:nrow(Tree$edge)){
      for (j in 1:ncol(Tree$edge)){
        matricePos[i,j]=dist.nodes(Tree)[Tree$edge[i,j],nbSpTot+1]
      }
    }
    ListNumBranche <- list()
    while (length(ListNumBranche) <= 1){
      N <- runif(1, min = matricePos[branche,1], max= matricePos[branche,2])
      for (i in 1:nrow(matricePos)){
        if((matricePos[i,1] < N) & (matricePos[i,2] > N)){
          ListNumBranche[length(ListNumBranche)+1] <- i
        }
      }
    }
    #Il faut qu'il y ai au moins deux point de coupure à un temps t pour réaliser un déplacement (un acctepeur et un donneur)
    
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
        b$edge.length=c(dist.nodes(Tree)[noeudOut,nbSpTot+1]-N,0)
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
        c$edge.length=c(dist.nodes(Tree)[noeudOut,nbSpTot+1]-N,0)
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
      tree3=Tree
      tree3$tip.label=Names
      treebinded = bind.tree(tree3, bindedBranch, where=noeudIns, position = dist.nodes(Tree)[noeudIns,nbSpTot+1]-N)
      treebinded = drop.tip(treebinded,"out")
      treebinded = drop.tip(treebinded,"ancienNoeud")
      return(treebinded)
}

