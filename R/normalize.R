#' normalize
#' 
#' This function normalizes the 2WR matrix (or any 2D matrix) according to the
#' species (rows) or to the genes (columns).
#' 
#' @param mat A matrix
#' @param what Character string indicating whether the matrix should be
#' normalized and how. If what="none", the matrix is not normalized (the
#' default), if what="species", the matrix is normalized so that the difference
#' between species is increased, and if what="genes", the matrix is normalized
#' so that the difference between genes is increased. Normalization consists
#' in dividing either each row or each columns by its median.
#' @return A normalized matrix
#' @export
normalize <- function(mat, what = "none") {
  if (what == "species") mat <- apply(mat, 2, function(x) {x / median(x)})
  else if (what == "genes") mat <- t(apply(mat, 1, function(x) {x / median(x)}))
  else if (what == "none") {
    mat <- mat
    #cat("The matric is unchanged.\n")
  } else print ("WARNING! Error in the kind of scaling you want! No scaling applied.")
  return(mat)
}

