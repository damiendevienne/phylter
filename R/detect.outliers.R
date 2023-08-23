# Detection of all types of outliers in the 2WR matrix  
# Three functions here.


#' Detection of outliers in 1D and 2D data
#' 
#' Functions to detect outliers, in matrices or in arrays. 
#' 
#' These functions detect outliers either in matrices or in arrays. For the method 
#' to be adapted to skewed data, as is the case here, the outlier detection 
#' method used is the adjusted Tukey proposed by Hubert and Vandervieren (2008).
#' 
#' @describeIn detect.outliers detect outliers in 2D matrix
#' 
#' @param X 2D matrix (gene x species) obtained with the \code{Dist2WR()} function.
#' @param k strength of outlier detection. High values (typically 3) 
#' correspond to less outliers detected than with lower ones (e.g. 1.5).
#' @param test.island should islands of outliers be treated as such. 
#' If TRUE (the default), only the highest value in an island of outliers is
#' removed. This prevents erroneous outliers due to hicthiking effect to be removed.
#' @param normalizeby Should the 2D matrix be normalized prior to outlier detection, and how.
#' Can be "row" (the default),"col" or "none". Normalization is done by dividing columns or rows 
#' by their median.
#' @return \code{detect.outliers}: A matrix with outliers detected in the 2D matric. Each row \code{x} contains the 
#' gene (\code{x[1]}) where the species (\code{x[2]}) is an outlier.  
#' @references de Vienne D.M., Ollier S. et Aguileta G. (2012) Phylo-MCOA: 
#' A Fast and Efficient Method to Detect Outlier Genes and Species
#' in Phylogenomics Using Multiple Co-inertia Analysis. Molecular 
#' Biology and Evolution 29 : 1587-1598.
#' @references Hubert, M. and Vandervieren, E. (2008). An adjusted boxplot for skewed
#' distributions. Computational Statistics and Data Analysis, 52, 5186-5201. 
#' @examples
#' # Get the initial gene x species matrix
#' # from the carnivora dataset
#' data(carnivora) 
#' mat <- phylter(carnivora, InitialOnly = TRUE, parallel = FALSE)$WR
#' 
#' # detect outliers in this matrix
#' outliers<-detect.outliers(mat)
#' outliers$cells # matrix where each row represents one cell in the input matrix
#' 
#' @importFrom stats dist quantile IQR
#' @export
detect.outliers <- function(X, k = 3, test.island=TRUE, normalizeby="row") {
  # heatmap(X, scale="none",Rowv=NA, Colv=NA)
  # scan()
  MAT <- X
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
  #normalized matrix so that each column (gene) is divided by its median
  if (normalizeby=="row") MATspgn <- normalize(X,"genes")
  else {
    if (normalizeby=="col") MATspgn <- normalize(X,"species")
    else MATspgn <- X
  }
  tabgn.TF<-outl.adj.tuckey(MATspgn,k)
  testFALSE<-tabgn.TF+0 
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
      RESULT$cells<-cbind(match(res[,2][-1], colnames(X)),match(res[,1][-1], rownames(X)))
    }
    else { ##we return all cells viewed as outliers (more violent...)
      cells<-which(testFALSE!=0,arr.ind = TRUE, useNames=FALSE)
      RESULT$cells<-cells[,c(2,1), drop=FALSE]
    }
  }
  else RESULT$cells<-NULL
#   plot(cbind(rep(1:ncol(X), each=nrow(X)), rep(1:nrow(X))), cex=array(X))
#   points(RESULT$cells, col="red")
# #  title(nrow(RESULT$cells))
#   scan() 
  RESULT$outgn<-NULL
  RESULT$outsp<-NULL
  return(RESULT)
}

#' @describeIn detect.outliers 
#' detects outliers in 1D array 
#' @param arr Array of values, typically the weight of each gene matrix (alpha values).
#' @param nbsp Number of species in the analysis
#' @return \code{detect.outliers.array}: An array listing the outliers detected (if any)  
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
