#' Filter phylogenomics datasets
#' 
#' Detection and filtering out of outliers in a list of trees 
#' or a list of distance matrices.
#' 
#' @param X A list of phylogenetic trees (phylo object) or a list 
#' of distance matrices. Trees can have different number of leaves and matrices
#' can have different dimensions. If this is the case, missing values are imputed. 
#' @param bvalue If X is a list of trees, nodes with a support below 'bvalue' will be collapsed
#' prior to the outlier detection.
#' @param distance If X is a list of trees, type of distance used to compute the 
#' pairwise matrices for each tree. Can be "patristic" (sum of branch lengths separating tips, the default)
#' or nodal (number of nodes separating tips). The "nodal" option should only be used if all species 
#' are present in all genes.
#' @param k Strength of outlier detection. The higher this value the less outliers
#' detected (see details).
#' @param k2 Same as k for complete gene outlier detection. To preserve complete genes from 
#' being discarded, k2 can be increased. By default, k2 = k. (see above) 
#' @param Norm Should the matrices be normalized prior to the complete analysis and how. If "median" (the default), matrices are divided by their median, if 
#' "mean" they are divided by their mean, if "none", no normalization if performed. Normalizing ensures that fast-evolving 
#' (and slow-evolving) genes are not treated as outliers. Normalization by median is a better choice as it is less sensitive to outlier values. 
#' @param Norm.cutoff Value of the median (if \code{Norm="median"}) or the mean (if
#' \code{Norm="mean"}) of phylogenetic distance matrices below which genes are simply discarded from the analysis. This
#' prevents dividing by 0, and allows getting rid of genes that contain mostly branches
#' of length 0 and are therefore uninformative anyway. Discarded genes, if any, are listed in 
#' the output ({out$DiscardedGenes}).
#' @param gene.names List of gene names used to rename elements in X. If NULL (the default), 0
#' elements are named 1,2,..,length(X). 
#' @param test.island If TRUE (the default), only the highest value in
#' an 'island' of outliers is considered an outlier. This prevents non-outliers hitchhiked by outliers
#' to be considered outliers themselves. 
#' @param verbose If TRUE (the default), messages are written during the filtering process to get information
#' of what is happening
#' @param stop.criteria The optimisation stops when the gain in concordance between matrices between round \code{n} and round \code{n+1} is smaller
#' than this value. Default to 1e-5.
#' @param InitialOnly Logical. If TRUE, only the Initial state of the data is computed. 
#' @param normalizeby Should the gene x species matrix be normalized prior to outlier detection, and how.
#' @return A list of class 'phylter' with the 'Initial' (before filtering) and 'Final' (after filtering) states, 
#' or a list of class 'phylterinitial' only, if InitialOnly=TRUE. The function also returns the list of DiscardedGenes, if any. 
#' @examples
#' data(carnivora)
#'
#' # using default paramaters
#' res<-phylter(carnivora) #perform the phylter analysis
#' res # brief summary of the analysis
#' res$DiscardedGenes # list of genes discarded prior to the analysis
#' res$Initial #See all elements prior to the analysis
#' res$Final #See all elements at the end of the analysis
#' res$Final$Outliers #Print all outliers detected
#' 
#' 
#' # Change the call to phylter  to use nodal distances, instead of patristic: 
#' res<-phylter(carnivora, distance="nodal")
#' 
#' @importFrom utils tail combn
#' @importFrom stats hclust as.dist median
#' @importFrom graphics plot
#' @export

phylter<-function(X, bvalue=0, distance="patristic", k=3, k2=k, Norm="median", Norm.cutoff=1e-3, gene.names=NULL, test.island=TRUE, verbose=TRUE, stop.criteria=1e-5, InitialOnly=FALSE, normalizeby="row") {
	ReplaceValueWithCompromise<-function(allmat, what, compro, lambda) {
		for (i in 1:length(allmat)) {
			whatsp<-what[what[,1]==i,2]
			if (length(whatsp)>0) {
				allmat[[i]][whatsp,]<-compro[whatsp,]*lambda[i]
				allmat[[i]][,whatsp]<-compro[,whatsp]*lambda[i]
			}
		}
		return(allmat)
	}
	ReoderBackTheCells<-function(DetectedCells, OrderUsed) {
		DetectedCells$cells[,2]<-OrderUsed[DetectedCells$cells[,2]]
		if (!is.null(DetectedCells$outsp)) DetectedCells$outsp<-OrderUsed[DetectedCells$outsp]
		return(DetectedCells)
	}
	FindNewCells<-function(DetectedCells, KnownCells) {
		NewCells<-NULL
		if (is.null(KnownCells)&!is.null(DetectedCells)) NewCells<-DetectedCells
		if (!is.null(KnownCells)&is.null(DetectedCells)) NewCells<-NULL
		if (!is.null(KnownCells)&!is.null(DetectedCells)) {
			#we paste the DetectedCells at the end of KnownCells, reusing KnownCells to save memory
			KnownCells<-rbind(KnownCells, DetectedCells)
			test<-tail(duplicated(KnownCells), nrow(DetectedCells))
			NewCells<-DetectedCells[test==FALSE,,drop=FALSE]
			if (nrow(NewCells)==0) NewCells<-NULL
		}
		return(NewCells)
	}
	CompareBeforeAfter<-function(InitialMatrices, AllOutliers, sp.order, which="all") {
		if (is.null(AllOutliers)) Out<-NULL
		else {
			if (which=="all") {
				cell.exists<-apply(AllOutliers, 1, function(x, sp) is.element(sp.order[x[2]],colnames(InitialMatrices[[x[1]]])), sp=sp.order)
				AllOutliers<-AllOutliers[cell.exists,]
				Out<-cbind(names(InitialMatrices[AllOutliers[,1]]), sp.order[AllOutliers[,2]])
			}
			if (which=="complete") {
				Out<-NULL
				speciesintables<-table(unlist(lapply(InitialMatrices, rownames)))
				speciesinoutliers<-table(AllOutliers[,2])
				Out$ComplOutSP<-names(which(speciesinoutliers==speciesintables[match(names(speciesinoutliers), names(speciesintables))]))
				genesintables<-unlist(lapply(InitialMatrices, ncol))
				genesinoutliers<-table(AllOutliers[,1])		
				Out$ComplOutGN<-names(which(genesinoutliers==genesintables[match(names(genesinoutliers), names(genesintables))]))
			}
		}
		return(Out)
	}

	#NEW
	X.clean<-PreparePhylterData(X, bvalue, distance, Norm, Norm.cutoff, gene.names, verbose)
	Xsave<-X.clean$Xsave
	matrices<-X.clean$matrices
	discardedgenes<-X.clean$discardedgenes #list pof discarded genes
	discardedmatrix<-X.clean$discardedmatrix #matrix of discarded cells (same format as Outliers at the end)
	#END NEW
	RES<-DistatisFast(matrices)
	WR<-Dist2WR(RES)

	Initial<-NULL
	Initial$mat.data<-Xsave
	Initial$WR<-WR
	Initial$RV<-RES$RVmat
	Initial$weights<-RES$alpha
	Initial$compromise<-RES$compromise
	Initial$F<-RES$F
	##New
	Initial$matrices<-matrices
	Initial$PartialF<-RES$PartialF

	if (InitialOnly) {
		class(Initial)<-c("phylterinitial", "list")
		return(Initial)
	}

	maxWR<-max(WR) ##for plotting purpose
	if (verbose) cat(paste("Initial score: ",round(RES$quality, digits=ceiling(abs(log10(stop.criteria)))), "", sep=""))
	VAL<-RES$quality #First Quality value
	CELLSREMOVED<-NULL
	continue<-TRUE
	lastLoop<-FALSE #when switching to TRUE, we do the last loop where we search for complete gene outliers.
	##We integrate this in the loop of cell outliers to keep possible the fact of not removing complete outliers detectected if they do not improve the quality of the compromise (very unlikely??)
	##OPTIMIZATION
	while(continue) {
		# cluster rows to be able to detect 'islands'
		OrderWRrow<-hclust(as.dist(RES$compromise))$order

		# reorder rows according to row clustering
		WR.reordered<-WR[OrderWRrow,]		
		if (!lastLoop) {
			# detect cell outliers
			CellsToRemove<-detect.outliers(WR.reordered, k=k,test.island=test.island, normalizeby=normalizeby)
			# reorder to previous order
			CellsToRemove<-ReoderBackTheCells(CellsToRemove, OrderWRrow)
		}
		else {
			CellsToRemove<-detect.outliers.array(RES$alpha, nrow(RES$compromise), k=k2)
		}
		NewCellsToRemove<-FindNewCells(CellsToRemove$cells, CELLSREMOVED)
		if (verbose & !is.null(nrow(NewCellsToRemove))) {
			cat(paste("\n   ",nrow(NewCellsToRemove)," new cells to remove "))
		}
		if (!is.null(NewCellsToRemove)) {
			CELLSREMOVED.new<-rbind(CELLSREMOVED, NewCellsToRemove) 
			matrices.new<-ReplaceValueWithCompromise(matrices, CELLSREMOVED.new, RES$compromise, RES$lambda)
			RES.new<-DistatisFast(matrices.new)
			VAL.new<-c(VAL, RES.new$quality)
			if (verbose) cat(paste("-> New score: ",round(RES.new$quality, digits=ceiling(abs(log10(stop.criteria)))), sep=""))
#			if (verbose) plot(VAL, type="o")				
			gain<-VAL.new[length(VAL.new)]-VAL.new[length(VAL.new)-1]
#			if (gain<0) {
			if (gain<stop.criteria) {
				if (verbose) cat(" -> NO")
				if (verbose) cat("\n => Gain too small")
				if (lastLoop) {
					continue<-FALSE #c'est vraiment la fin
					if (verbose) cat("  ->  Stopping optimization.")

				}
				else {
					if (verbose) cat ("  ->  Checking for complete gene outliers")
					lastLoop<-TRUE ##
				}
				
			}
			else { #on continue. 
				if (verbose) cat(" -> OK")
				RES<-RES.new
				matrices<-matrices.new
				CELLSREMOVED<-CELLSREMOVED.new
				VAL<-VAL.new
				WR<-Dist2WR(RES)
				#si lastLoop était vrai, il redevient faux car on se demande si on peut réaméliorer maintenant améliorer.
				if (lastLoop) lastLoop<-FALSE
			}

		}
		else {
			if (verbose) cat ("\n => No more outliers detected")
			if (lastLoop) {
				cat ("  ->  STOPPING OPTIMIZATION")
				break ##usefull?
				continue<-FALSE #c'est vraiment la fin
			}
			else {
				if (verbose) cat ("  ->  Checking for complete gene outliers")
				lastLoop<-TRUE

			}
		}
	}

	# Prepare return object
	Final<-NULL
	Final$WR<-WR
	Final$RV<-RES$RVmat
	Final$weights<-RES$alpha
	Final$compromise<-RES$compromise
	Final$F<-RES$F
	Final$PartialF<-RES$PartialF
	Final$species.order<-colnames(RES$compromise)
	Final$AllOptiScores<-VAL

	# ### NEW - 10/03/2022
	# ### WE MUST MODIFY CELLSREMOVED TO REINTEGRATE THE GENES
	# ### THAT WERE REMOVED BECAUSE OF MEDIAN or MEAN=0 ISSUE. INDEED,
	# ### if gene n is removed, gene n+1 becomes n, n+2 becomles n+1, etc.
	# ### We need to reorganise these cells to put the correct identifiers.
	# ### we do that as follows:
	# print(CELLSREMOVED)
	# if (length(ZeroLengthMatrices)>0) {
	# 	for (i in 1:length(ZeroLengthMatrices)) {
	# 		CELLSREMOVED[CELLSREMOVED[,1]>=ZeroLengthMatrices[i],1]<-CELLSREMOVED[CELLSREMOVED[,1]>=ZeroLengthMatrices[i],1]+1
	# 	}
	# 	#and then we add as outliers these genes that had been removed
	# 	#and all there species
	# 	CELLSREMOVED<-rbind(CELLSREMOVED, cbind(rep(ZeroLengthMatrices, each=nrow(WR)),rep(1:nrow(WR),length(ZeroLengthMatrices))))

	# }
	# ### END OF THIS NEW PART 

	Final$CELLSREMOVED<-CELLSREMOVED 
	#MAYBE DO NOT RETURN THIS!?	
	Final$Outliers<-CompareBeforeAfter(Xsave, CELLSREMOVED, Final$species.order, which="all")
	Final$CompleteOutliers<-CompareBeforeAfter(Xsave, Final$Outliers, Final$species.order, which="complete")
	#store the way the function was called
	##New 
	Final$matrices<-matrices
	Final$Discarded<-discardedmatrix #organized like Outliers, but for genes discarded at the very beginning
	
	call<-list(call=match.call(), bvalue=bvalue, distance=distance, k=k, k2=k2, Norm=Norm, Norm.cutoff=Norm.cutoff, gene.names=gene.names, test.island=test.island, verbose=verbose, stop.criteria=stop.criteria, InitialOnly=InitialOnly, normalizeby=normalizeby)

	Result<-list(Initial=Initial, Final=Final, DiscardedGenes=discardedgenes, call=call)
	class(Result)<-c("phylter", "list")
	class(Result$Initial)<-c("phylterinitial", "list")
	class(Result$Final)<-c("phylterfinal", "list")

	if (verbose) {
		cat("\n--------\n")
		print(summary(Result))
	}
	return(Result)
}
