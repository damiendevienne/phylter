# New implementation of Distatis, much faster when only conserving a few axes.

#' Fast implementation the multivariate analysis method Distatis
#' 
#' New implementation of the DISTATIS method for K matrices of dimension IxI.
#' This version of Distatis is faster than the original one because only the minimum required number
#' of eigenvalues and eigenvectors is calculated. The difference in speed 
#' is particularly visible when the number of matrices is large.  
#' 
#' @param matrices A list of K distance matrices, all of the same dimension (IxI).
#' @param factorskept Number of factors to keep for the computation of the factor 
#' @param parallel Should the matrix products computations be parallelized? Default to TRUE. 
#' scores of the observations.
#' @return Returns a list containing:
#' \itemize{
#'  \item \code{F}: projected coordinates of species in the compromise
#'  \item \code{PartialF}: list of projected coordinates of species for each gene.
#' 	\item \code{alpha}: array of length K of the weight associated to each matrix.
#'  \item \code{lambda}: array of length K of the normalization factors used for each matrix. lambda=1 always.
#' 	\item \code{RVmat}: a KxK matrix with RV correlation coefficient computed between all
#'  pairs of matrices.
#' 	\item \code{compromise}: an IxI matrix representing the best compromise between all matrices. This matrix is the weighted average of all K matrices, using 'alpha' as
#' a weighting array.
#' 	\item \code{quality}: the quality of the compromise. This value is between 0 and 1 
#' 	\item \code{matrices.dblcent}: matrices after double centering
#' describes how much of the variance of the K matrices is captured by the compromise.   
#' }
#' @references Abdi, H., Valentin, D., O'Toole, A.J., & Edelman, B. (2005).
#' DISTATIS: The analysis of multiple distance matrices.
#' Proceedings of the IEEE Computer Society: International
#' Conference on Computer Vision and Pattern Recognition_.  (San
#' Diego, CA, USA). pp. 42-47.
#' @examples
#' # Get a list of matrices 
#' # from the carnivora dataset
#' data(carnivora) 
#' matrices<-phylter(carnivora, InitialOnly=TRUE)$matrices
#' 
#' # Perform a Distatis analysis on these matrices: 
#' distatis<-DistatisFast(matrices)
#' 
#' #distatis is a list with multiple elements: 
#' distatis$alpha #weigh of each matrix (how much it correlates with others)
#' distatis$RVmat #RV matrix: correlation of each matrix with each other
#' distatis$compromise # distance matrix with "average" pairwise distance between species in matrices
#' # etc.
#' 
#' @importFrom RSpectra eigs_sym
#' @importFrom Rfast Crossprod 
#' @export
DistatisFast<-function(matrices, factorskept=2, parallel=TRUE) {
	GetCmat <- function(OrderedMatrices, RV = TRUE, parallel) {
	    CP2.diag <-do.call(cbind, lapply(OrderedMatrices, diag))
	    CP2.upper <- do.call(cbind, lapply(OrderedMatrices, function(x) x[upper.tri(x)]))
#	    C <- ifelse(parallel, Crossprod(CP2.diag,CP2.diag) + 2 * Crossprod(CP2.upper,CP2.upper), crossprod(CP2.diag,CP2.diag) + 2 * crossprod(CP2.upper,CP2.upper))
	    C <- Crossprod(CP2.diag,CP2.diag) + 2 * Crossprod(CP2.upper,CP2.upper)	       
	    if (RV) {
	        laNorm = sqrt(2 * colSums(CP2.upper^2) + colSums(CP2.diag^2))
	        C = C/outer(laNorm, laNorm)
	    }
	    rownames(C) <- colnames(C) <- names(OrderedMatrices)
	    return(C)
	} 
	DblCenterDist <- function(Y) {
		##ACCELERATION POSSIBLE ? VOIR BICENTER DANS ADE4
		# nI = nrow(Y)
		# CentMat = diag(nI) - (1/nI) * matrix(1, nI, nI)
		# S = -(1/2) * (CentMat %*% Y %*% CentMat)
		# dimnames(S)<-dimnames(Y)
		# return(S)	 
		col.mean<-colMeans(Y)
		row.mean<-rowMeans(Y)
    	Y <- sweep(Y, 2, col.mean)
	    Y <- sweep(Y, 1, row.mean)
    	Y <- Y + mean(col.mean)
    	return(-Y/2)

	}
	MFAnormCP <- function(Y) {
	    e1 = eigs_sym(Y, 1, opts=list(retvec=FALSE), which="LA")$values
	    Ynormed = Y/e1
	    return(Ynormed)
	}
	GetLambdaForNorm <- function(Y) {
	    e1 = eigs_sym(Y, 1, opts=list(retvec=FALSE), which="LA")$values
	    return(e1)
	}
	Sp<-rownames(matrices[[1]])
	Gn<-names(matrices)
	nbSp<-length(Sp)
	nbGn<-length(Gn)
	### Compute weights
	matrices.dblcent<-lapply(matrices, DblCenterDist)
	# if (Norm) matrices.dblcent<-lapply(matrices.dblcent, MFAnormCP) ##normalize is asked
	lambda<-rep(1,nbGn)
	RVmat<-GetCmat(matrices.dblcent, parallel=parallel)
	FirstEigenVector<-eigs_sym(RVmat, 1, which = "LM")
	alpha <- FirstEigenVector$vectors[, 1]/sum(FirstEigenVector$vectors[, 1])
	quality<-FirstEigenVector$values/nbGn
	### Compute compromise matrix (C) and its projection (Splus)
	WeightedMatrices<-sapply(1:nbGn, function(x,MAT,weight) MAT[[x]]*weight[x],MAT=matrices.dblcent, weight=alpha, simplify=FALSE)
	Splus<-Reduce('+',WeightedMatrices)
	# compromise<-Reduce('+',WeightedMatrices.initial)
	dimnames(Splus)<-list(Sp,Sp)
	s<-diag(Splus)
	#### dé-centrage ####
	compromise<-sweep(sweep(-2*(Splus),1,s,"+"),2,s,"+")
	## ca fois me mambda c'est OK.
	### Keep few (=factorskept) axes and project individual matrices
	Nom2Factors<-paste("Factor", 1:factorskept)
	eigenSplus = eigs_sym(Splus, factorskept) ##the 2 is the number of dim we really keep.
	eigenSplus$SingularValues<-sqrt(abs(eigenSplus$values))
	F<-t(apply(eigenSplus$vectors, 1, "*", t(t(eigenSplus$SingularValues))))
	rownames(F) <- Sp
	colnames(F) <- Nom2Factors
	Proj <- t(apply(eigenSplus$vectors, 1, "*", 1/t(t(eigenSplus$SingularValues))))
	colnames(Proj) <-  paste("Factor", 1:ncol(F))
	rownames(Proj) <- Sp
	PartialF = lapply(matrices.dblcent, function(x,y) x %*% y, y=Proj)
	return(list(F=F, PartialF=PartialF, alpha=alpha, lambda=lambda, RVmat=RVmat, compromise=compromise, quality=quality, matrices.dblcent=matrices.dblcent))
}

