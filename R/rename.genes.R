# This function renames objects in a list (typcally trees or matrices). 
# If no name is given, names are given as number between 1 and the number of elements in X



#' rename.genes
#' 
#' This function names or renames trees in a list of trees.
#' 
#' 
#' @param X A list of trees or matrices
#' @param gene.names List of names to assign to the elements of X. Must be of the same length as length(X). 
#' If NULL (the default) the object are numbered 1,2,...,length(X).
#' @return X with name assigned to each element.
rename.genes <- function(X, gene.names = NULL) {
  if(!is.null(gene.names)) {
    names(X) <- gene.names
  } else {
    if(is.null(names(X))) {
      names(X) <- as.character(c(1:length(X)))
    }
  }
  return(X)
}
