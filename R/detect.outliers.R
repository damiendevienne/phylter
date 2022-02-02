# Detection of all types of outliers in the 2WR matrix  
# Three functions here.


#' Detect outliers
#' 
#' Functions to detect outliers, either complete or cell (see details). 
#' 
#' These functions detect outliers in the 2WR matrix. The method 
#' used is adapted to skewed data, as is the case here. Prior to outlier
#' detection, the 2WR matrix is normalized by row and by column
#' and the outliers are detected from the complete matrix using the 
#' adjusted Tukey proposed by Hubert and Vandervieren (2008): 
#' Hubert, M. and Vandervieren, E. (2008). An adjusted boxplot for skewed
#' distributions. Computational Statistics and Data Analysis, 52, 5186-5201. 
#' The different types of outliers (cell and complete) are described in:
#' de Vienne D.M., Ollier S. et Aguileta G. (2012) Phylo-MCOA: 
#' A Fast and Efficient Method to Detect Outlier Genes and Species
#' in Phylogenomics Using Multiple Co-inertia Analysis. Molecular 
#' Biology and Evolution 29 : 1587-1598.
#' 
#' @describeIn detect.outliers a simple wrapper
#' to call the two other functions described, and return complete outliers
#' if any, or cell outliers if no complete outliers exist.
#' 
#' @param mat2WR the 2WR matrix obtained with the Dist2WR function.
#' @param k the strength of outlier detection. High values (typically 3) 
#' correspond to less outliers detected than with lower ones (e.g. 1.5).
#' @param test.island should islands of outliers be treated as such. 
#' Default to FALSE. If TRUE, only the highest value in the island is 
#' removed (This option may be deprecated soon). 
#' @param old Should the old detection method be used instead (default to FALSE).
#' @return A list of outliers.
#' @references de Vienne D.M., Ollier S. et Aguileta G. (2012) Phylo-MCOA: 
#' A Fast and Efficient Method to Detect Outlier Genes and Species
#' in Phylogenomics Using Multiple Co-inertia Analysis. Molecular 
#' Biology and Evolution 29 : 1587-1598.
#' Hubert, M. and Vandervieren, E. (2008). An adjusted boxplot for skewed
#' distributions. Computational Statistics and Data Analysis, 52, 5186-5201. 
#' @importFrom mrfDepth medcouple
#' @importFrom stats dist quantile IQR
#' @export
detect.outliers <- function(mat2WR, k = 3, test.island=TRUE, old=FALSE) {
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
  #NEW
  outl.adj.tuckey<-function(x,k) {
    mc<-medcouple(array(x))
    return(x > (quantile(x)[4] + k*exp(3*mc)*IQR(x) + 1e-10))
  }
  #END NEW
  outl.sub <- function(x, k) {
    return(x > quantile(x)[4] + k * IQR(x) + 1e-10)
  }
  if (old) {
    MATspgn <- normalize(mat2WR, "genes") * normalize(mat2WR, "species")
    testspgn1 <- apply(MATspgn, 2, outl.sub, k = k)
    testspgn2 <- t(apply(MATspgn, 1, outl.sub, k = k))
    testspgn <- testspgn1 * testspgn2 
    testFALSE <- testspgn
  }
  else {
    #normalized matrix so that each column (gene) is divided by its median
    MATspgn <- normalize(mat2WR,"species")
    tabgn.TF<-outl.adj.tuckey(MATspgn,k)
    testFALSE<-tabgn.TF+0 
  }
  #
  RESULT<-NULL
  #
  if (sum(testFALSE) > 0) {
    if (test.island==TRUE) {
      out.list <- apply(testFALSE, 2, detect.island)
      genes <- colnames(testFALSE)
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
      RESULT$cells<-cells[,c(2,1), drop=FALSE]
    }
  }
  else RESULT$cells<-NULL
#   plot(cbind(rep(1:ncol(mat2WR), each=nrow(mat2WR)), rep(1:nrow(mat2WR))), cex=array(mat2WR))
#   points(RESULT$cells, col="red")
# #  title(nrow(RESULT$cells))
#   scan() 
  RESULT$outgn<-NULL
  RESULT$outsp<-NULL
  return(RESULT)
}

#' @describeIn detect.outliers detects
#' if array (here the correlations between gene matrices) contains outlier. This
#' is the last step of the phylter process.
#' @param arr Array of values, typically the weight of each gene matrix (alpha values).
#' @param nbsp Number of species in the analysis
#' @export
detect.outliers.array <- function(arr, nbsp, k = 3) {
  RES<-NULL 
  #new outlier detection method:
  outl.adj.tuckey<-function(x,k) {
    mc<-medcouple(array(x))
    return(x > (quantile(x)[4] + k*exp(3*mc)*IQR(x) + 1e-10))
  }
  out.genes<-which(outl.adj.tuckey(-arr, k=k)) #we put negative(arr) because we want the **smallest** outliers
  if (length(out.genes)>0) {
    RES$outgn<-out.genes
    cells.outgn<-cbind(rep(out.genes, each=nbsp),rep(1:nbsp, length(out.genes))) 
    RES$cells<-cells.outgn
    RES$outsp<-NULL
  }
  else {
    RES$outgn<-NULL
    RES$outsp<-NULL
    RES$cells<-NULL
  }
  return(RES)
}
