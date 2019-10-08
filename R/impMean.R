# Impute missing datas in matrices by computing the mean of non-empty corresponding entries

#' impMean
#' 
#' Impute missing data in a list of matrices. Matrices are first given the same dimension, 
#' then missing entries are computed by taking the average value in non-missing corresponding
#' entries in all matrices.
#' 
#' 
#' @param matrices A list of distance matrices. 
#' @return Returns a list of matrices with same dimensions, with rows and columns 
#' in the same order and missing data (if any) imputed.
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
  qual <- 0
  sp.per.mat<-lapply(matrices, rownames)
  listsp <- unique(unlist(sp.per.mat))
  length.indiv <- unlist(lapply(sp.per.mat, length))
  tmp <- rep(1,length(length.indiv))
  if (sum(length(listsp)-length.indiv)>0) qual<-1 
  if (qual == 1) {
    matrices.extended<-lapply(matrices, AddColAndRow, allrowsandcol=listsp) #grows matrices and add NA to missing cells
    MEANMAT<-Reduce('+',lapply(matrices.extended, function(x) replace(x,is.na(x),0)))/Reduce("+", lapply(matrices.extended, Negate(is.na)))
    ALL<-lapply(matrices.extended, ReplaceMissingValue, matrix2=MEANMAT)
  }
  else { ##it is a simple reordering of the columnss and rows
    ALL<-lapply(matrices, ReorderColAndRow, allrowsandcol=listsp)
  }
  names(ALL) = names(matrices)
  return(ALL)
}

