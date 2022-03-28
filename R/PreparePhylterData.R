#' Prepare data for phylter analysis
#' 
#' Prepare datasets for the \code{phylter} function. Detects possible issues, 
#' discads genes if necessary, imputes missing data if any, and reorders row- and col-names.
#' For internal usage mostly.
#'
#' @param X A list of phylogenetic trees (phylo object) or a list 
#' of distance matrices. Trees can have different number of leaves and matrices
#' can have different dimensions. If this is the case, missing values are imputed. 
#' @param bvalue If X is a list of trees, nodes with a support below 'bvalue' will be collapsed
#' prior to the outlier detection.
#' @param distance If X is a list of trees, type of distance used to compute the 
#' pairwise matrices for each tree. Can be "patristic" (sum of branch lengths separating tips, the default)
#' or nodal (number of nodes separating tips).
#' @param Norm Should the matrices be normalized and how. If "median", matrices are divided by their median, if 
#' "mean" they are divided by their mean, if "none", no normalization if performed. Normalizing ensures that fast-evolving 
#' (and slow-evolving) genes are not treated as outliers. Normalization by median is less sensitive to outlier values
#' but can lead to errors if some matrices have a median value of 0. 
#' are not considered outliers. 
#' @param Norm.cutoff Value of the median (if Norm="median") or the mean (if
#' Norm="mean") below which matrices are simply discarded from the analysis. This
#' prevents dividing by 0, and getting rid of genes that contain mostly branches
#' of length 0 and are therefore uninformative anyway. 
#' @param gene.names List of gene names used to rename elements in X. If NULL (the default), 
#' elements are named 1,2,..,length(X). 
#' @param verbose If TRUE (the default), messages are written during the filtering process to get information
#' of what is happening
#' @return A list of class 'phylter' with the 'Initial' (before filtering) and 'Final' (after filtering) states, 
#' or a list of class 'phylterinitial' only, if InitialOnly=TRUE. 
#' @examples 
#' data(carnivora)
#' # transform trees to a named list of matrices with same dimensions
#' # and identical row and column orders and names
#' carnivora.clean<-PreparePhylterData(carnivora)
#' 
#'  
#' @importFrom utils tail combn
#' @importFrom stats hclust as.dist median
#' @importFrom graphics plot
#' @export

PreparePhylterData<-function(X, bvalue=0, distance="patristic", Norm="median",Norm.cutoff=1e-6,gene.names=NULL, verbose=TRUE) {

	WarningMessage.GenesDiscarded<-function(nb,nam) {
		if (nb>1) {
			cat("\n########################################\n")
			cat("###           WARNING!               ###\n")
			cat(paste("### ",nb, " genes were removed prior",paste(rep(" ",8-nchar(nb)), collapse=""),"###\n",sep=""))
			cat("### to the analysis because they had ###\n")
			cat("### too many zero-length branches    ###\n")
			cat("########################################\n")
			cat(paste("Dicarded genes: ",paste(paste("\"",nam, "\"", sep=""),collapse=", "),"\n", sep=""))

		}
		else {
			cat("\n########################################\n")
			cat("###           WARNING!               ###\n")
			cat(paste("### ",nb, " gene was removed prior",paste(rep(" ",10-nchar(nb)), collapse=""),"###\n",sep=""))
			cat("### to the analysis because it had   ###\n")
			cat("### too many zero-length branches    ###\n")
			cat("########################################\n")
			cat(paste("Dicarded genes: ",paste(paste("\"",nam, "\"", sep=""),collapse=", "),"\n", sep=""))
		}
	}
 

	if (is.null(names(X))) X<-rename.genes(X, gene.names=gene.names)
	if (inherits(X[[1]], "phylo")) matrices <- trees2matrices(X, distance = distance, bvalue = bvalue)
	else matrices<-X
	Xsave<-matrices #Xsave contains the original matrices
	matrices <- impMean(matrices) ##impute missing values with mean (if any). This also sorts rows and columns, thus this step cannot be removed.
	if (verbose) {
		cat(paste("\nNumber of Genes:    ", length(matrices), "\n", sep=""))
		cat(paste("Number of Species:  ", nrow(matrices[[1]]), "\n", sep=""))
		cat("--------\n")
	}
	## NEW (25/01/2022). Doing this, genes that havehigh values globally are not treated as outliers, 
	## but long branches remain outliers (median is not too affected by outliers).
	## NEW March 10 2022 
	##remove matrices that have mean or median close to 0
	ZeroLengthMatrices<-NULL
	discardedmatrix<-NULL
	if (Norm %in% c("median","mean")) {
		if (Norm=="median") {
			AllMe<-unlist(lapply(matrices, median))
			matrices<-lapply(matrices, function(x) x/median(x))
			} else if (Norm=="mean") {
				AllMe<-unlist(lapply(matrices, mean))
				matrices<-lapply(matrices, function(x) x/mean(x))
			}	
		#remove matrices with median or mean < Norm.cutoff
		ZeroLengthMatrices<-which(AllMe<=Norm.cutoff)
		if (length(ZeroLengthMatrices)>0) {
			#for later use, get list of all discarded species for all discarded genes
			SpeciesPerDiscardedGenes<-lapply(Xsave[ZeroLengthMatrices], rownames)
			matrices<-matrices[-ZeroLengthMatrices]
			#and remove those matrices from Xsave as well
			Xsave<-Xsave[-ZeroLengthMatrices]
			if (verbose) {
					WarningMessage.GenesDiscarded(length(ZeroLengthMatrices), names(ZeroLengthMatrices))
				}
			###prepare the discarded matrix (col1=genes, col2=species)
			discardedmatrix<-cbind(rep(names(ZeroLengthMatrices),unlist(lapply(SpeciesPerDiscardedGenes,length))), unname(unlist(SpeciesPerDiscardedGenes)))
		}
	}	
	discardedgenes<-names(ZeroLengthMatrices)

	return(list(matrices=matrices, Xsave=Xsave, discardedgenes=discardedgenes, discardedmatrix=discardedmatrix))
}

