#' normalize
#' 
#' This function normalizes the 2WR matrix (or any matrix) according to the
#' species (rows) or to the genes (columns). This normalization is described
#' in de Vienne M.D., Ollier S. et Aguileta G. (2012) Phylo-MCOA: 
#' A Fast and Efficient Method to Detect Outlier Genes and Species
#' in Phylogenomics Using Multiple Co-inertia Analysis. Molecular Biology 
#' and Evolution 29 : 1587 – 1598)
#' 
#' @param mat A matrix
#' @param what Character string indicating whether the matrix should be
#' normalized and how. If what="none", the matrix is not normalized (the
#' default), if what="species", the matrix is normalized so that the difference
#' between species is increased, and if what="genes", the matrix is normalized
#' so that the difference between genes is increased. Normalization consists
#' in dividing each entry by either the mean of each row or of each column.
#' @return A normalized matrix
#' @export
normalize <- function(mat, what = "none") {
  if (what == "species") mat <- apply(mat, 2, function(x) {x / mean(x)})
  else if (what == "genes") mat <- t(apply(mat, 1, function(x) {x / mean(x)}))
  else if (what == "none") {
    mat <- mat
    #cat("The matric is unchanged.\n")
  } else print ("WARNING! Error in the kind of scaling you want! No scaling applied.")
  return(mat)
}
