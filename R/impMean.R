# Impute missing datas in matrices by computing the mean of non-empty corresponding entries

#' Imputation of missing values in a collection of matrices
#' 
#' Impute missing data in a list of matrices. Matrices are first given the same dimension, 
#' then missing entries are computed by taking the average value in non-missing corresponding
#' entries in all matrices.
#' 
#' 
#' @param matrices A list of distance matrices. 
#' @return Returns a list of matrices with same dimensions, with rows and columns 
#' in the same order and missing data (if any) imputed.
#' @examples 
#' data(carnivora)
#' matrices<-phylter(carnivora, InitialOnly=TRUE)$matrices
#' 
#' # remove n species randomly (n between 1 and 5) in each matrix (to mimic missing data)
#' fun<-function(mat) {
#'  species2remove<-sample(1:nrow(mat),sample(1:5,1))
#'  mat<-mat[-species2remove,-species2remove]
#'  return(mat)
#' }
#' matrices.missing<-lapply(matrices, fun)
#' #check that all matrices have now different dimensions: 
#' lapply(matrices.missing, dim)
#' # Impute data to get back to the same dimensions
#' matrices.ok<-impMean(matrices.missing)
#' lapply(matrices.ok, dim) #all dimensions are now identical. Missing data have been imputed. 
#'
#' @importFrom reshape2 melt
#' @export
impMean <- function(matrices) {
  AddColAndRow<-function(matrix, allrowsandcol) {
    newmat<-matrix(nrow=length(allrowsandcol), ncol=length(allrowsandcol), dimnames=list(allrowsandcol,allrowsandcol))
    newmat[rownames(matrix),colnames(matrix)]<-matrix
    newmat
  }
  ReorderColAndRow<-function(matrix, allrowsandcol) {
    newmat<-matrix[allrowsandcol,allrowsandcol]
    newmat
  }
  ReplaceMissingValue<-function(matrix1,matrix2) { ##this function will replace NA in matrix1 by the corresponding value in matrix2
    matrix1[is.na(matrix1)]<-matrix2[is.na(matrix1)]
    matrix1
  }
  ImputeFromClosestNeighbors<-function(sp1,sp2,mat) {
    closest2sp1<-sp1
    closest2sp2<-sp2
    i<-1
    j<-1
    while(is.na(mat[sp1,closest2sp2])) {
      i<-i+1
      closest2sp2<-names(sort(mat[sp2,])[i])
    }
    while(is.na(mat[sp2,closest2sp1])) {
      j<-j+1
      closest2sp1<-names(sort(mat[sp1,])[j])
    }
    imputeddist<-mean(mat[sp1,closest2sp2],mat[sp2,closest2sp1]) 
    imputeddist
  }

  qual <- 0
  sp.per.mat<-lapply(matrices, rownames)
  listsp <- unique(unlist(sp.per.mat))
  length.indiv <- unlist(lapply(sp.per.mat, length))
  tmp <- rep(1,length(length.indiv))
  if (sum(length(listsp)-length.indiv)>0) qual<-1 
  if (qual == 1) {
    matrices.extended<-lapply(matrices, AddColAndRow, allrowsandcol=listsp) #grows matrices and add NA to missing cells
    MEANMAT<-Reduce('+',lapply(matrices.extended, function(x) replace(x,is.na(x),0)))/Reduce("+", lapply(matrices.extended, Negate(is.na)))
    #If NA are present in this matrix, it means that some species were never found together in any tree. If so we need to find a way to imputethe distance by looking at the distance of their close neighbors (and print a warning). 
    if(anyNA(MEANMAT)) {
      tmpM<-melt(MEANMAT)
      pairs2correct<-tmpM[is.na(tmpM$value),1:2]
      newval<-apply(pairs2correct,1, function(x,y) ImputeFromClosestNeighbors(x[1],x[2],y), y=MEANMAT)
      MEANMAT[as.numeric(rownames(pairs2correct))]<-newval
      cat("\nWARNING: some pairs of species are never found together in trees. \nMean distances are thus infered using their closest neighbors.\n")
    }
    ALL<-lapply(matrices.extended, ReplaceMissingValue, matrix2=MEANMAT)
  }
  else { ##it is a simple reordering of the columnss and rows
    ALL<-lapply(matrices, ReorderColAndRow, allrowsandcol=listsp)
  }
  names(ALL) = names(matrices)
  return(ALL)
}

