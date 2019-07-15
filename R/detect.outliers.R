# Detection of all types of outliers in the 2WR matrix  
# Three functions here.


#' Detect outliers
#' 
#' Functions to detect outliers, either complete or cell (see details). 
#' 
#' These methods are similar to those described in phylo-MCOA for detecting
#' outliers. These types of outliers are described in:
#' de Vienne M.D., Ollier S. et Aguileta G. (2012) Phylo-MCOA: 
#' A Fast and Efficient Method to Detect Outlier Genes and Species
#' in Phylogenomics Using Multiple Co-inertia Analysis. Molecular 
#' Biology and Evolution 29 : 1587-1598.
#' 
#' @describeIn detect.outliers a simple wrapper
#' to call the two other functions described, and return complete ouliers
#' if any, or cell outliers if no complete outliers exist.
#' 
#' @param mat2WR the 2WR matrix obtained with the Dist2WR function.
#' @param k the strength of outlier detection. High values (typically 3) 
#' correspond to less outliers detected than with lower ones (e.g. 1.5).
#' @param thres threshold above which genes or species are considered as
#' complete outliers. 0.3 means that a gene (resp. species) is a
#' complete outlier if it is detected as outlier for more than 30\% of 
#' its species (resp. genes).
#' @param test.island should islands of outliers be treated as such. 
#' Default to TRUE, i.e. only the highest value in the island is 
#' removed. See details. 
#' @param keep.species should species be prevented from being removed. This
#' is useful if all species are important in the dataset. Default to TRUE
#' (all species are kept).
#' @param outlier.detection.method Method used to detect outliers from the 2WR matrix. Default to 1.
#' @return A list of outliers.
#' @export
detect.outliers<-function(mat2WR, k=3, thres=0.3, test.island=TRUE, keep.species=TRUE, outlier.detection.method=1) {
  #Test for complete outliers
  CELLS<-NULL
  if (outlier.detection.method==1) CompOutl <- detect.complete.outliers(mat2WR, k = k, thres = thres, keep.species=keep.species)
  if (outlier.detection.method==2) CompOutl <- detect.complete.outliers2(mat2WR, k = k, thres = thres, keep.species=keep.species)
  if (nrow(CompOutl$cells)>0) {
    CELLS$outgn<-CompOutl$outgn
    CELLS$outsp<-CompOutl$outsp
    CELLS$cells<-CompOutl$cells
  }
  else {
  if (outlier.detection.method==1) CellOutl <- detect.cell.outliers1(mat2WR, k = k, test.island=test.island)
  if (outlier.detection.method==2) CellOutl <- detect.cell.outliers2(mat2WR, k = k, test.island=test.island)
    CELLS$outgn<-NULL
    CELLS$outsp<-NULL
    CELLS$cells<-CellOutl$cells
  }
  return(CELLS)
}

#' @describeIn detect.outliers detects
#' if complete outliers exist in the 2WR matrix. See details.
#' @export
detect.complete.outliers <- function(mat2WR, k = 3, thres = 0.3, keep.species=TRUE) {
  RES<-NULL
  outl.sub <- function(x, k) {
    return(x > quantile(x)[4] + k * IQR(x) + 1e-10)
    ##note: the 1e-10 ,is because when all values are similar except one, the first one is considered as equal to the third quartile... May be a bug in quantile function?
  }
  tabgn <- normalize(mat2WR, "genes")
  tabgn.TF <- t(apply(tabgn, 1, outl.sub, k = k))
  tabgn.TF<-tabgn.TF+0
  score.genes <- apply(tabgn.TF, 2, mean)
  out.genes <- which(score.genes > thres)
#  RES$outgn <- names(out.genes)
  RES$outgn <- out.genes
  cells.outgn<-cbind(rep(out.genes, each=nrow(mat2WR)),rep(1:nrow(mat2WR), length(out.genes))) 
  if (!keep.species) {
    tabsp <- normalize(mat2WR, "species")
    tabsp.TF <- apply(tabsp,2, outl.sub, k = k)
    tabsp.TF<-tabsp.TF+0
    score.species <- apply(tabsp.TF, 1, mean)
    out.species <- which(score.species > thres)
    # RES$outsp <- names(out.species)
    RES$outsp <- out.species
    cells.outsp<-cbind(rep(1:ncol(mat2WR), length(out.species)), rep(out.species, each=ncol(mat2WR)))
  }
  else {
    RES$outsp<-NULL 
    cells.outsp<-NULL
  }
  RES$cells <- rbind(cells.outgn, cells.outsp)  #Col 1 = gn; col 2 = sp
  dimnames(RES$cells)<-NULL
  return(RES)
}

#' @describeIn detect.outliers detects
#' if complete outliers exist in the 2WR matrix. New version.
#' @export
detect.complete.outliers2 <- function(mat2WR, k = 3, thres = 0.3, keep.species=TRUE) {
  RES<-NULL
  outl.sub <- function(x, k) {
    return(x > quantile(x)[4] + k * IQR(x) + 1e-10)
    ##note: the 1e-10 ,is because when all values are similar except one, the first one is considered as equal to the third quartile... May be a bug in quantile function?
  }

  mat2WR_normalized<-t(apply(mat2WR,1,function(x) x/mean(x))) #nroamlize by row
  testspgn<-apply(mat2WR_normalized,2,function(x) outl.sub(x,k=k)) + 0
  score.genes <- apply(testspgn, 2, mean)
  out.genes <- which(score.genes > thres)
  RES$outgn <- out.genes
  cells.outgn<-cbind(rep(out.genes, each=nrow(mat2WR)),rep(1:nrow(mat2WR), length(out.genes))) 
  if (!keep.species) {
    score.species <- apply(testspgn, 1, mean)
    out.species <- which(score.species > thres)
    # RES$outsp <- names(out.species)
    RES$outsp <- out.species
    cells.outsp<-cbind(rep(1:ncol(mat2WR), length(out.species)), rep(out.species, each=ncol(mat2WR)))
  }
  else {
    RES$outsp<-NULL 
    cells.outsp<-NULL
  }
  RES$cells <- rbind(cells.outgn, cells.outsp)  #Col 1 = gn; col 2 = sp
  dimnames(RES$cells)<-NULL
  return(RES)
}


#' @describeIn detect.outliers detects
#' if cell outliers exist in the 2WR matrix.
#' @importFrom stats dist quantile IQR
#' @importFrom utils combn
#' @export
detect.cell.outliers <- function(mat2WR, k = 3, test.island=TRUE) {
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
  #
  RESULT<-NULL
  #
  if (sum(testFALSE) > 0) {
    if (test.island==TRUE) {
      out.list <- apply(testspgn, 2, detect.island)
      genes <- colnames(testspgn)
      res <- c(NA,NA)
      for (i in 1:length(out.list)) {
        if (!is.null(out.list[[i]])) {
          for (j in 1:length(out.list[[i]])) { ##for each "island"
            if (length(out.list[[i]][[j]]) == 1) res <- rbind(res, c(out.list[[i]][[j]], genes[i]))
            if (length(out.list[[i]][[j]]) > 1) {
              vals <- MATspgn[out.list[[i]][[j]], genes[i]]
              multi = c(names(vals)[vals == max(vals)])
              if(length(multi) > 1){
                x = cbind(multi, rep(genes[i],length(multi)))
                res <- rbind(res, x)
              }
              else{
                res <- rbind(res, c(multi, genes[i]))
              }
            }
          }
        }
      }
      colnames(res) <- c("Species", "Genes")
      RESULT$cells<-cbind(match(res[,2][-1], colnames(mat2WR)),match(res[,1][-1], rownames(mat2WR)))
    }
    else { ##we return all cells viewed as outliers (more violent...)
      cells<-which(testFALSE!=0,arr.ind = TRUE, useNames=FALSE)
      RESULT$cells<-cells[,c(2,1)]
    }
  }
  return(RESULT)
}


#' @describeIn detect.outliers detects
#' if cell outliers exist in the 2WR matrix. New version.
#' @importFrom stats dist quantile IQR
#' @importFrom utils combn
#' @export
detect.cell.outliers2 <- function(mat2WR, k = 3, test.island=TRUE) {
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
  ### IN THIS NEW VERSION, WE NORMALIZE FIRST SO THAT ALL SPECIES HAVE THE SAME MEAN
  ### THEN WE DETECT THE OUTLIERS BY COLUMNS, INDEPENDENTLY FOR EACH. 
  ### 

  mat2WR_normalized<-t(apply(mat2WR,1,function(x) x/mean(x))) #nroamlize by row
  testspgn<-apply(mat2WR_normalized,2,function(x) outl.sub(x,k=k)) + 0
  #
  RESULT<-NULL
  #
  if (sum(testspgn) > 0) {
    if (test.island==TRUE) {
      out.list <- apply(testspgn, 2, detect.island)
      genes <- colnames(testspgn)
      res <- c(NA,NA)
      for (i in 1:length(out.list)) {
        if (!is.null(out.list[[i]])) {
          for (j in 1:length(out.list[[i]])) { ##for each "island"
            if (length(out.list[[i]][[j]]) == 1) res <- rbind(res, c(out.list[[i]][[j]], genes[i]))
            if (length(out.list[[i]][[j]]) > 1) {
              vals <- mat2WR[out.list[[i]][[j]], genes[i]]
              multi = c(names(vals)[vals == max(vals)])
              if(length(multi) > 1){
                x = cbind(multi, rep(genes[i],length(multi)))
                res <- rbind(res, x)
              }
              else{
                res <- rbind(res, c(multi, genes[i]))
              }
            }
          }
        }
      }
      colnames(res) <- c("Species", "Genes")
      RESULT$cells<-cbind(match(res[,2][-1], colnames(mat2WR)),match(res[,1][-1], rownames(mat2WR)))
    }
    else { ##we return all cells viewed as outliers (more violent...)
      cells<-which(testspgn!=0,arr.ind = TRUE, useNames=FALSE)
      RESULT$cells<-cells[,c(2,1)]
    }
  }
  return(RESULT)
}
