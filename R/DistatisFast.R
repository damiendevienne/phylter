# New implementation of Distatis, much faster when only conserving a few axes.

#' DistatisFast
#' 
#' New implementation of the DISTATIS method for K matrices of dimension IxI.
#' This version of distatis is faster than the original one because only the minimum required number
#' of eigenvalues and eigenvectors is calculated. The difference in speed 
#' is particularly visible when the number of matrices is large.  
#' 
#' @param matrices A list of K distance matrices, all of the same dimension (IxI).
#' @param Norm Should the matrices be normalized. If TRUE (the default), 
#' each matrix is normalized such that its first eigenvalie is equal to one.
#' @param factorskept Number of factors to keep for the computation of the factor 
#' scores of the observations.
#' @return Returns a list contanining 
#' \itemize{
#' 	\item 'alpha': array of length K of the weight associated to each matrix.
#'  \item 'lambda': array of length K of the normalization factors used for each matrix 
#'  if Norm=TRUE.
#' 	\item 'RVmat': a KxK matrix with RV correlation coefficient computed between all
#'  pairs of matrices.
#' 	\item 'compromise': an IxI matrix representing the best compromise between all 
#' matrices. This matrix is the weighted average of all K matrices, using 'alpha' as
#' a weighting array.
#' 	\item 'quality': the quality of the compromise. This value between 0 and 1 
#' describes how much of the variance of the K matrices is captured by the compromise.   
#' }
#' @references Abdi, H., Valentin, D., O'Toole, A.J., & Edelman, B. (2005).
#' DISTATIS: The analysis of multiple distance matrices.
#' Proceedings of the IEEE Computer Society: International
#' Conference on Computer Vision and Pattern Recognition_.  (San
#' Diego, CA, USA). pp. 42-47.
#' @importFrom RSpectra eigs_sym
#' @export
DistatisFast<-function(matrices, Norm=TRUE, factorskept=2) {
	GetCmat <- function(OrderedMatrices, RV = TRUE) {
		CP2<-do.call(cbind, lapply(OrderedMatrices, array))
		C <- crossprod(CP2) #faster than using t(CP2) %*% CP2
	    if (RV) {
	        laNorm = sqrt(apply(CP2^2, 2, sum))
	        C = C/(t(t(laNorm)) %*% laNorm)
	    }
	    rownames(C) <- colnames(C) <- names(OrderedMatrices)
	    return(C)
	}
	GetCmat2 <- function(OrderedMatrices, RV = TRUE) {
		dim1<-nrow(OrderedMatrices[[1]])
		index<-unlist(sapply(1:dim1,function(x, n) (x:n)+(x-1)*n, n=dim1))		
		CP2<-do.call(cbind, lapply(OrderedMatrices, function(x,y) array(x)[y], y=index))
		C <- crossprod(CP2) #faster than using t(CP2) %*% CP2
	    if (RV) {
	        laNorm = sqrt(apply(CP2^2, 2, sum))
	        C = C/(t(t(laNorm)) %*% laNorm)
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
	if (Norm) {
		lambda<-lapply(matrices.dblcent, GetLambdaForNorm)
		matrices.dblcent<-Map("/", matrices.dblcent, lambda)
		lambda<-unname(unlist(lambda))
	}
	else {
		lambda<-rep(1,nbGn)
	}
	RVmat<-GetCmat(matrices.dblcent)
	FirstEigenVector<-eigs_sym(RVmat, 1, which = "LM")
	alpha <- FirstEigenVector$vectors[, 1]/sum(FirstEigenVector$vectors[, 1])
	quality<-FirstEigenVector$values/nbGn
	### Compute compromise matrix (C) and its projection (Splus)
	WeightedMatrices<-sapply(1:nbGn, function(x,MAT,weight) MAT[[x]]*weight[x],MAT=matrices.dblcent, weight=alpha, simplify=FALSE)
	WeightedMatrices.initial<-sapply(1:nbGn, function(x,MAT,weight) MAT[[x]]*weight[x],MAT=matrices, weight=alpha, simplify=FALSE)

	Splus<-Reduce('+',WeightedMatrices)
	# compromise<-Reduce('+',WeightedMatrices.initial)
	dimnames(Splus)<-list(Sp,Sp)
	s<-diag(Splus)
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

	return(list(F=F, PartialF=PartialF, alpha=alpha, lambda=lambda, RVmat=RVmat, compromise=compromise, quality=quality))
}
