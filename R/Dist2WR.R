
# \code{Dist2WR} computes the 2WR matrix from the results obtaind with \code{DistatisFast}
# (the fast version of distatis)



#' Dist2WR
#' 
#' \code{Dist2WR} computes the 2WR matrix from the results obtaind with \code{DistatisFast}
# (the fast version of distatis).
#' 
#' 
#' @param Distatis is the output of the fonction \code{DistatisFast}.
#' @return The 2WR matrix: a gene x species matrix. Each cell represents a 
#' gene/species pair, whose value represents the distance between (i) the position of this species
#' in this gene tree and (ii) the average position of this species in all other gene trees.
#' 
Dist2WR <- function(Distatis) {
  F<-Distatis$F
  PartialF<-Distatis$PartialF
  DISTS<-lapply(PartialF, function(f,fbar) apply((f-fbar)^2,1,function(x) sqrt(sum(x))), fbar=F)
  matrixWR2<-do.call(cbind, DISTS)
  return(matrixWR2)
}
