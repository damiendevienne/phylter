# trees2matrices changes a list of trees into a list of matrices



#' trees2matrices
#' 
#' trees2matrices changes a list of trees into a list of matrices.
#' 
#' 
#' @param trees A list of gene trees in multiphylo format.
#' @param distance A method to generate distance matrices. It could be "nodal"
#' to establish that the distance between two species is the number of nodes
#' that separate them. Or "patristic" (default) if the distance between two
#' species is be the sum of branch lengths between them.
#' @param bvalue This argument is only used if trees contain bootstrap values.
#' It determines under what bootstrap values the nodes should be collapsed.
#' Value 0 (the default) means that no nodes are collapsed.
#' @return return a list of distance matrices
#' @importFrom ape Nnode Ntip di2multi compute.brlen
#' @importFrom stats cophenetic
#' @export
trees2matrices <- function(trees, distance = "patristic", bvalue = 0) {
  correction <- function(mat){
    mat<-mat-1
    diag(mat)<-0
    return(mat)
  }
  list.trees <- list()
  for (i in 1:length(trees)) {
    tree <- trees[[i]]
    if (distance == "nodal") {
      if (bvalue != 0) {
        if (!is.null(tree$node.label)) {
          l <- 1:Nnode(tree)
          indices.nodes <- l[as.numeric(tree$node.label) < bvalue] + Ntip(tree)
          if (length(indices.nodes) > 0) {
            for (j in 1:length(indices.nodes)) {
              tree$edge.length[tree$edge[,1] == indices.nodes[j]] <- 1e-10
            }
          }
          tree <- di2multi(tree, tol = 1e-9)
        }
        else {
          #tree <- di2multi(tree, tol = bvalue)
          cat ("-------- WARNING!! There are no bootstraps values in your trees or your bvalue parameter does not correspond to the values in the trees ! -------")
        }
      }
      tree.brlen <- compute.brlen(tree, 1)
    }
    else if (distance == "patristic") {
      if (bvalue != 0) {
        if (!is.null(tree$node.label)) {
          l <- 1:Nnode(tree)
          indices.nodes <- l[as.numeric(tree$node.label) < bvalue] + Ntip(tree)
          if (length(indices.nodes) > 0) {
            for (j in 1:length(indices.nodes)) {
              tree$edge.length[tree$edge[,1] == indices.nodes[j]] <- 1e-10
            }
          }
          tree <- di2multi(tree, tol = 1e-9)
        }
        else {
          #tree <- di2multi(tree, tol = bvalue)
          cat ("-------- WARNING!! There are no bootstraps values in the trees ! -------")
        }
      }
      tree.brlen <- tree
    }
    list.trees[[i]] <- tree.brlen
  }
  TRS <- lapply(list.trees, cophenetic)
  if (distance == "nodal"){
    TRS <- lapply(TRS, correction)
  }
  if (!is.null(names(trees))) {
    names(TRS) <- names(trees)
  }
  return(TRS)
}
