#' Compute gene x species matrix from the result of Distatis
#' 
#' \code{Dist2WR} computes the 2WR matrix from the results obtaind with \code{DistatisFast}
#' (the fast version of distatis). For internal use mostly.
#' 
#' 
#' @param Distatis output of the fonction \code{DistatisFast}.
#' @return A matrix (gene=rows x species=col). Each cell represents a 
#' gene/species pair, whose value represents the distance between (i) the position of this species
#' in this gene tree and (ii) the average position of this species in all other gene trees.
#' @examples
#' data(carnivora)
#' matrices<-phylter(carnivora, InitialOnly=TRUE, parallel=FALSE)$matrices
#' ds<-DistatisFast(matrices, parallel = FALSE)
#' WR<-Dist2WR(ds) #returns the gene x species matrix
#' 
#' @export
Dist2WR <- function(Distatis) {
  F<-Distatis$F
  PartialF<-Distatis$PartialF
  DISTS<-lapply(PartialF, function(f,fbar) apply((f-fbar)^2,1,function(x) sqrt(sum(x))), fbar=F)
  matrixWR2<-do.call(cbind, DISTS)
  return(matrixWR2)
}
