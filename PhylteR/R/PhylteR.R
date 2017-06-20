# trees2matrices changes a list of trees into a list of matrices

trees2matrices <- function(trees, distance = "patristic", bvalue = 0) {
  correction <- function(mat){
    for (i in 1:nrow(mat)){
      for (j in 1:ncol(mat)){
        if (i != j) {mat[i,j] <- mat[i,j] - 1}
      }
    }
    return(mat)
  }
  list.trees<-list()
  for (i in 1:length(trees)) {
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
  TRS<-lapply(list.trees, cophenetic)
  if (distance == "nodal"){
    TRS <- lapply(TRS,correction)
  }
  if (!is.null(names(trees))) {
    names(TRS) <- names(trees)
  }
  return(TRS)
}

# this function permits the user to add gene names (as list) to the trees. If no names is given, genes are numeroted from 1 to the number of genes.

rename.genes <- function(trees, gene.names = NULL) {
  if(!is.null(gene.names)) {
    names(trees) <- gene.names
  } else {
    if(is.null(names(trees))) {
      names(trees) <- as.character(c(1:length(trees)))
    }
  }
  return(trees)
}

# mat2Dist applies distatis on a list of distance matrices

mat2Dist <- function(matrices, Norm = "NONE") {
  # transform the list of matrices to a cube
  row <- rownames(matrices[[1]])
  for (i in 1:length(matrices)) {
    matrices[[i]] <- matrices[[i]][row, row]
  }
  genesNumber <- length(matrices)
  speciesNumber <- nrow(matrices[[genesNumber]])
  TheVeryBigCube <- array(0, c(speciesNumber, speciesNumber, genesNumber))
  for (i in 1:genesNumber) {
    TheVeryBigCube[, , i] <- matrices[[i]]
  }
  rownames(TheVeryBigCube) <- rownames(matrices[[1]])
  colnames(TheVeryBigCube) <- colnames(matrices[[1]])
  # Apply distatis on the cube and keep genes names in distatis results
  Distatis <- distatis(TheVeryBigCube, Norm = Norm)
  dimnames(Distatis$res4Splus$PartialF)[[3]] <- names(matrices)
  return(Distatis)
}

# Imputing missing data in matrices with missMDA package

impPCA.multi <- function(matrices, ncp = 3, center = FALSE, scale = FALSE, maxiter = 10000) {
  geneNames <- list()
  # Create a matrix with every species to fill : GrandeMatrice
  species<-unique(unlist(lapply(matrices, rownames)))
  nbsp<-length(species) #ca évite de le recalculer plein de fois.
  grandeMatrice=matrix(nrow = length(species), ncol = nbsp)
  rownames(grandeMatrice) = species
  colnames(grandeMatrice) = species
  # Create a list of matrices names with missing data
  dimMat<-unlist(lapply(matrices, nrow))
  geneNames<-names(which(dimMat<nbsp))
  matrices2 = matrices
  if (length(geneNames) != 0) {
    # For each matrix with missing data, we fill the GrandeMatrice
    for (i in 1:length(geneNames)) {
      row <- rownames(matrices2[[geneNames[[i]]]])
      grandeMatrice2 <- grandeMatrice
      grandeMatrice2[row, row] <- matrices2[[geneNames[[i]]]][row,row]
      matrices2[[geneNames[[i]]]] <- grandeMatrice2
    }
    matrices3 <- matrices2
    row <- rownames(matrices2[[1]])
    for (i in 1:length(matrices2)) {
      matrices3[[i]] <- as.vector(matrices2[[i]][row, row][upper.tri(matrices2[[i]][row, row])])
    }
    mat <- do.call(cbind,matrices3)
    # estimating missing data
    matIPCA <- imputePCA2(mat, center = center, scale = scale, maxiter = maxiter)
    matIPCA <- matIPCA$completeObs
    matricesFT <- list()
    for (i in 1:length(matrices2)) {
      matrices[[i]] <- matrix(nrow = length(row), ncol = length(row))
      rownames(matrices[[i]]) <- row
      colnames(matrices[[i]]) <- row
      data=as.vector(matIPCA[, i])
      matrices[[i]][upper.tri(matrices2[[i]])] <- data
      matricesFT[[i]] <- t(matrices[[i]])
      matrices[[i]][lower.tri(matrices2[[i]])] <- matricesFT[[i]][lower.tri(matrices2[[i]])]
      for (j in 1:ncol(matrices[[i]])) {
        matrices[[i]][j, j] <- 0
      }
    }
  }
  return(matrices)
}

# imputePCA2 function from missMDA package but with the possibility to not center the data. And with no negative values possibly imputed.

imputePCA2 <- function (X, ncp = 2, center = FALSE, scale = FALSE, method = c("Regularized", "EM"), row.w = NULL, coeff.ridge = 1, threshold = 1e-6, seed = NULL,nb.init = 1, maxiter = 1000) {
  impute <- function (X, ncp = 4, center = FALSE, scale = FALSE, method = NULL, threshold = 1e-6,seed = NULL, init = 1, maxiter = 1000, row.w = NULL, coeff.ridge = 1) {
    moy.p <- function(V, poids) {
      res <- sum(V * poids, na.rm = TRUE) / sum(poids[!is.na(V)])
    }
    ec <- function(V, poids) {
      res <- sqrt(sum(V ^ 2 * poids,na.rm = TRUE) / sum(poids[!is.na(V)]))
    }
    nb.iter <- 1
    old <- Inf
    objective <- 0
    if (!is.null(seed)) {
      set.seed(seed)
      }
    X <- as.matrix(X)
    ncp <- min(ncp, ncol(X), nrow(X) - 1)
    missing <- which(is.na(X))
    mean.p <- apply(X, 2, moy.p, row.w)
    if (center != TRUE){
      mean.p <- rep(0, ncol(X))
    }
    Xhat <- t(t(X) - mean.p)
    et <- apply(Xhat, 2, ec, row.w)
    if (scale) Xhat <- t(t(Xhat) / et)
    if (any(is.na(X))) Xhat[missing] <- 0
    if (init > 1) Xhat[missing] <- rnorm(length(missing)) # random initialization
    fittedX <- Xhat
    if (ncp == 0) nb.iter = 0
    while (nb.iter > 0) {
      Xhat[missing] <- fittedX[missing]
      if (scale) Xhat = t(t(Xhat) * et)
      Xhat <- t(t(Xhat) + mean.p)
      if (center == TRUE){
        mean.p <- apply(Xhat, 2, moy.p, row.w)
      }
      Xhat <- t(t(Xhat) - mean.p)
      et <- apply(Xhat, 2, ec, row.w)
      if (scale) Xhat <- t(t(Xhat) / et)
      svd.res <- svd.triplet(Xhat, row.w = row.w, ncp = ncp)
      sigma2  <- nrow(X) * ncol(X) / min(ncol(X), nrow(X)-1) * sum((svd.res$vs[-c(1:ncp)] ^ 2) / ((nrow(X) - 1) * ncol(X) - (nrow(X) - 1) * ncp - ncol(X) * ncp + ncp ^ 2))
      sigma2 <- min(sigma2 * coeff.ridge, svd.res$vs[ncp + 1] ^ 2)
      if (method == "em") sigma2 <-0
      lambda.shrinked = (svd.res$vs[1:ncp] ^ 2 - sigma2) / svd.res$vs[1:ncp]
      fittedX = tcrossprod(t(t(svd.res$U[, 1:ncp, drop = FALSE] * row.w) * lambda.shrinked), svd.res$V[, 1:ncp, drop = FALSE])
      fittedX <- fittedX / row.w
      # No negative value
      fittedX[fittedX < 0] = 0
      diff <- Xhat-fittedX
      diff[missing] <- 0
      objective <- sum(diff ^ 2 * row.w)
      criterion <- abs(1 - objective / old)
      old <- objective
      nb.iter <- nb.iter + 1
      if (!is.nan(criterion)) {
        if ((criterion < threshold) && (nb.iter > 5)) nb.iter <- 0
        if ((objective < threshold) && (nb.iter > 5)) nb.iter <- 0
      }
      if (nb.iter > maxiter) {
        nb.iter <- 0
        warning(paste("Stopped after ", maxiter, " iterations"))
      }
    }
    if (scale) Xhat <- t(t(Xhat) * et)
    Xhat <- t(t(Xhat) + mean.p)
    completeObs <- X
    completeObs[missing] <- Xhat[missing]
    if (scale) fittedX <- t(t(fittedX) * et)
    fittedX <- t(t(fittedX) + mean.p)

    result <- list()
    result$completeObs <- completeObs
    result$fittedX <- fittedX
    return(result)
  }
  # Main program
  method <- match.arg(method, c("Regularized", "regularized", "EM", "em"), several.ok = T)[1]
  obj = Inf
  method <- tolower(method)
  if (ncp > min(nrow(X) - 2, ncol(X) - 1)) stop("ncp is too large")
  if (is.null(row.w)) row.w = rep(1, nrow(X)) / nrow(X)
  for (i in 1:nb.init){
    if (!any(is.na(X))) return(X)
    res.impute = impute(X, ncp = ncp, scale = scale,center = center, method = method, threshold = threshold, seed = if(!is.null(seed)){(seed * (i - 1))}else{NULL}, init = i, maxiter = maxiter, row.w = row.w, coeff.ridge = coeff.ridge)
    if (mean((res.impute$fittedX[!is.na(X)] - X[!is.na(X)]) ^ 2) < obj) {
      res <- res.impute
      obj <- mean((res.impute$fittedX[!is.na(X)] - X[!is.na(X)]) ^ 2)
    }
  }
  return(res)
}

# create 2WR matrix from distatis results

Dist2WR <- function(Distatis) {
  matrixWR2 <- matrix(nrow = dim(Distatis$res4Splus$PartialF)[[1]], ncol = dim(Distatis$res4Splus$PartialF)[[3]])
  colnames(matrixWR2) <- dimnames(Distatis$res4Splus$PartialF)[[3]]
  rownames(matrixWR2) <- dimnames(Distatis$res4Splus$PartialF)[[1]]

  for (i in 1:length(dimnames(Distatis$res4Splus$PartialF)[[3]])){
    for (j in 1:length(dimnames(Distatis$res4Splus$PartialF)[[1]])){
      x <- (Distatis$res4Splus$PartialF[dimnames(Distatis$res4Splus$PartialF)[[1]][j], , dimnames(Distatis$res4Splus$PartialF)[[3]][i]] - Distatis$res4Splus$F[dimnames(Distatis$res4Splus$PartialF)[[1]][j], ]) ^ 2
      matrixWR2[dimnames(Distatis$res4Splus$PartialF)[[1]][j], dimnames(Distatis$res4Splus$PartialF)[[3]][i]] <- sqrt(sum(x))
    }
  }
  return(matrixWR2)
}

# Fonction to plot 2WR matrix

 plot2WR <- function(matrixWR2) {
  WR <- normalize(matrixWR2)
  names <- list()
  names[[1]] <- "gene"
  names[[2]] <- "specie"
  names[[3]] <- "value"
  MAT <- matrix(nrow = length(WR), ncol = 3)
  colnames(MAT) <- names
  k <- 1
  for (i in 1:nrow(WR)) {
    for (j in 1:ncol(WR)) {
      MAT[k, 2] <- rownames(WR)[i]
      MAT[k, 1] <- colnames(WR)[j]
      MAT[k, 3] <- WR[i, j]
      k <- k + 1
    }
  }
  MAT <- as.data.frame(MAT)
  genes <- as.character(MAT$gene)
  species <- as.character(MAT$specie)
  distanceToRef <- as.numeric(as.character(MAT$value))

  pl <- ggplot(MAT, aes(genes, species, z = distanceToRef))
  pl <- pl + geom_tile(aes(fill = distanceToRef)) + theme_bw() + scale_fill_gradient(low = "white", high = "blue")
  pl <- pl + theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 13, color = "black"))
  pl <- pl + theme(axis.text.y = element_text(angle = 00, hjust = 1, size = 13, color = "black"))
  pl <- pl + theme(axis.title.y = element_text(size = rel(1.8), angle = 90))
  pl <- pl + theme(axis.title.x = element_text(size = rel(1.8), angle = 00))
  # pl + coord_fixed(ratio=1/5)
  return(pl)
}

# Suppress complete outiers (species or genes) in trees in order to detect cell outliers in a second time.
# Fonction from phylo_MCOA slitly modified to fit our method.

rm.gene.and.species <- function(trees, sp2rm, gn2rm) {
  gene.names <- list()
  for (i in 1:length(labels(trees))) {
    gene.names[i] <- labels(trees)[i]
  }
  if (length(sp2rm) > 0) {
    for (i in 1:length(trees)) {
      sp <- trees[[i]]$tip.label
      toremove <- intersect(sp, sp2rm)
      trees[[i]] <- drop.tip(trees[[i]], toremove)
    }
  }
  if (length(gn2rm) > 0) {
    genes2keep <- setdiff(gene.names, gn2rm)
    trees2 <- list()
    j <- 0
    for (i in 1:length(trees)) {
      if (is.element(gene.names[i], genes2keep)) {
        j <- j + 1
        trees2[[j]] <- trees[[i]]
      }
    }
    names(trees2) <- genes2keep
  } else {
    trees2 <- trees
  }
  return(trees2)
}

# Phylter Function to detect complete and cell outliers from a list of trees

PhylteR <- function(trees, distance = "patristic", k = 1.5, thres = 0.5, gene.names = NULL, Norm = "NONE", bvalue = 0) {
  if (is.list(trees)) {
    if (class(trees[[1]]) != "phylo") stop ("The trees should be in the \"phylo\" format!")
  }
  if (class(trees) == "character") {
    trees<-read.tree(trees)
  }
  if (!is.null(gene.names) && length(gene.names) != length(trees)) stop ("The number of gene names and the number of trees differ!")
  ##check for duplications
  check.dup<-lapply(trees, function(x) {x$tip.label[duplicated(x$tip.label)]})
  if (sum(unlist(lapply(check.dup, length)))>0) {
    cat ("-------- ATTENTION!! There are some duplicated species in some of your trees: -------")
    for (w in 1:length(trees)) {
      if (length(check.dup[[w]])>0) cat(paste("\n     - Species ", check.dup[[w]], " present more than once in tree ",w,"\n\n",sep=""))
    }
    stop ("Remove or rename duplicated species and try again.\n\n", call.=FALSE)
  }
  trees <- rename.genes(trees, gene.names = gene.names)
  RES <- NULL
  matrices <- trees2matrices(trees, distance = distance, bvalue = bvalue)
  matrices <- impPCA.multi(matrices)
  Dist <- mat2Dist(matrices, Norm = Norm)
  WR <- Dist2WR(Dist)
  CompOutl <- detect.complete.outliers(WR, k = k, thres = thres)
  if (length(CompOutl$outsp) > 0 || length(CompOutl$outgn) > 0) {
    TREESwithoutCompleteOutlierDist <- rm.gene.and.species(trees, CompOutl$outsp, CompOutl$outgn)
    matrices2 <- trees2matrices(TREESwithoutCompleteOutlierDist, distance = distance, bvalue = bvalue)
    matrices2 <- impPCA.multi(matrices2)
    Dist2 <- mat2Dist(matrices2, Norm = Norm)
    WR2 <- Dist2WR(Dist2)
    CellOutl2 <- detect.cell.outliers(WR2, k = k + 2)
    RES$Complete <- CompOutl
    RES$CellByCell <- CellOutl2
  } else {
    CellOutl2 <- detect.cell.outliers(WR,  k = k + 2)
    RES$Complete <- CompOutl
    RES$CellByCell <- CellOutl2
  }
  return(RES)
}

#This function normalizes the 2WR matrix (or any matrix) according to the species (rows) or to the genes (columns).

normalize <- function(mat, what = "none") {
  if (what == "species") mat <- apply(mat, 2, function(x) {x / mean(x)})
  else if (what == "genes") mat <- t(apply(mat, 1, function(x) {x / mean(x)}))
  else if (what == "none") {
    mat <- mat
    #cat("The matric is unchanged.\n")
  } else print ("WARNING! Error in the kind of scaling you want! No scaling applied.")
  return(mat)
}

# Detection of complete outliers from phylo-mcoa
detect.complete.outliers <- function(mat2WR, k = 1.5, thres = 0.5) {
  outl.sub <- function(x, k) {
    return(x > quantile(x)[4] + k * IQR(x) + 1e-10)
    ##note: the 1e-10 ,is because when all values are similar except one, the first one is considered as equal to the third quartile... May be a bug in quantile function?
  }
  tabgn <- normalize(mat2WR, "genes")
  tabgn.TF <- t(apply(tabgn, 1, outl.sub, k = k))
  tabsp <- normalize(mat2WR, "species")
  tabsp.TF <- apply(tabsp,2, outl.sub, k = k)
  tabgn.TF[tabgn.TF == FALSE] <- 0
  tabgn.TF[tabgn.TF == TRUE] <- 1
  tabsp.TF[tabsp.TF == FALSE] <- 0
  tabsp.TF[tabsp.TF == TRUE] <- 1
  score.genes <- apply(tabgn.TF, 2, function(x) {sum(x) / length(x)})
  score.species <- apply(tabsp.TF, 1, function(x) {sum(x) / length(x)})
  out.genes <- names(score.genes)[score.genes > thres]
  out.species <- names(score.species)[score.species > thres]
  RES <- NULL
  RES$mat2WR <- mat2WR
  #RES$thres <- thres
  #RES$allgn <- names(score.genes)
  #RES$allsp <- names(score.species)
  #RES$scoregn <- score.genes
  #RES$scoresp <- score.species
  #RES$TFgn <- score.genes>thres
  #RES$TFsp <- score.species>thres
  RES$outgn <- out.genes
  RES$outsp <- out.species
  return(RES)
}

# Function to detect cell outliers (species and genes)
detect.cell.outliers <- function(mat2WR, k = 3) {
  ans <- "y"
  if (ans == "y") {
    MAT <- mat2WR
    detect.island <- function(arr) {
      spi.names <- names(arr)
      spi <- 1:length(spi.names)
      names(spi) <- spi.names
      true.names <- names(arr)[arr == TRUE]
      if (length(true.names) == 1) {
        return(list(true.names))
      }
      else if (length(true.names) > 1) {
        true.i <- spi[true.names]
        res<-dist(true.i)
        table.i <- cbind(t(combn(attributes(res)$Labels, 2)), array(res))
        in.island<-NULL
        if (length(table.i[table.i[, 3] == "1", 3]) == 0) {
          in.island <- "nopair"
          list.i <- NULL
        }
        if (length(table.i[table.i[,3] == "1", 3]) == 1) {
          in.island <- table.i[table.i[,3] == "1", c(1, 2)]
          list.i <- list(in.island)
        }
        if (is.null(in.island)) {
          table.small <- table.i[table.i[,3] == "1",c(1,2)]
          list.i <- list()
          for (i in 1:nrow(table.small)) list.i[[i]] <- table.small[i, ]
          for (i in 1:(length(list.i) - 1)) {
            for (j in (i+1):length(list.i)) {
              if (length(intersect(list.i[[i]], list.i[[j]])) > 0) {
                list.i[[i]] <- c(list.i[[i]], list.i[[j]])
                list.i[[j]] <- "out"
                list.i[[i]] <- unique(list.i[[i]])
              }
            }
          }
          list.i2 <- list()
          w <- 0
          for (i in 1:length(list.i)) {
            if ((length(list.i[[i]]) > 1)&&(list.i[[i]][1] != "out")) {
              w <- w+1
              list.i2[[w]] <- list.i[[i]]
              in.island <- c(in.island, list.i[[i]])
            }
          }
          list.i <- list.i2
        }
        out.island <- as.list(setdiff(true.names, in.island))
        return(c(list.i, out.island))
      }
    }
    outl.sub <- function(x, k) {
      return(x > quantile(x)[4] + k * IQR(x) + 1e-10)
    }
    MATspgn <- normalize(mat2WR, "genes") * normalize(mat2WR, "species")
    testspgn1 <- apply(MATspgn, 2, outl.sub, k = k)
    testspgn2 <- t(apply(MATspgn, 1, outl.sub, k = k))
    testspgn <- testspgn1 * testspgn2
    testFALSE <- testspgn
    testFALSE[testFALSE == FALSE] <- 0
    testFALSE[testFALSE == TRUE] <- 1
    if (sum(testFALSE) > 0) {
      out.list <- apply(testspgn, 2, detect.island)
      genes <- colnames(testspgn)
      res <- c(NA,NA)
      for (i in 1:length(out.list)) {
        if (!is.null(out.list[[i]])) {
          for (j in 1:length(out.list[[i]])) {
            if (length(out.list[[i]][[j]]) == 1) res <- rbind(res, c(out.list[[i]][[j]], genes[i]))
            if (length(out.list[[i]][[j]]) > 1) {
              vals <- MATspgn[out.list[[i]][[j]], genes[i]]
              res <- rbind(res, c(names(vals)[vals == max(vals)], genes[i]))
            }
          }
        }
      }
      colnames(res) <- c("Species", "Genes")
      ##we construct the MATfinal
      MATfinal <- testspgn
      MATfinal[,] <- 0
      for (w in 2:nrow(res)) MATfinal[res[w, 1], res[w, 2]] <- 1
        RESULT <- NULL
        #RESULT$mat2WR <- mat2WR
        #RESULT$matspgn <- MATspgn
        #RESULT$matfinal <- MATfinal
        #RESULT$testFALSE <- testFALSE
        RESULT$outcell <- res[2:nrow(res), ]
      return(RESULT)
    } else return(NULL)
  } else cat("\n---Operation canceled by the user.---\n")
}
