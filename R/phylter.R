# Detect and filter outliers in a list of trees or distance matrices. 

#' filter phylogenomics datasets
#' 
#' Detection and filtering out of outliers in a list of trees 
#' or distance matrices.
#' 
#' @param X A list of phylogenetic trees (phylo object) or a list 
#' of distance matrices. Trees can have different number of leaves and matrices
#' can have different dimensions. If this is the case, missing values are imputed. 
#' @param bvalue If X is a list of trees, nodes with a support below 'bvalue' will be collapsed
#' prior to the outlier detection.
#' @param distance If X is a list of trees, type of distance used to compute the 
#' pairwise matrices for each tree. Can be "patristic" (sum of branch lengths separating tips, the default)
#' or nodal (number of nodes separating tips).
#' @param k Strength of outlier detection. The higher this value the less outliers
#' detected (see details).
#' @param k2 Same as k for complete gene outlier detection. To preserve complete genes from 
#' being discarded, k2 can be increased . By default, k2 = k. (see above) 
#' By default, k2=k. 
#' @param Norm Should the matrices be normalized. If TRUE (the default), 
#' each matrix is divided by its median. This ensures that fast-evolving genes
#' are not considered outliers. 
#' @param gene.names List of gene names used to rename elements in X. If NULL (the default), 
#' elements are named 1,2,..,length(X). 
#' @param test.island This should not be modified. If TRUE (the default), only the highest value in
#' an 'island' of outliers is considered an outlier. This prevents non-outliers hitchhiked by outliers
#' to be considered outliers themselves. 
#' @param verbose If TRUE (the default), messages are written during the filtering process to get information
#' of what is happening
#' @param stop.criteria The optimisation stops when the gain between round n and round n+1 is smaller
#' than this value. Default to 1e-5.
#' @param InitialOnly Logical. If TRUE, only the Initial state of teh data is computed. The optimization and 
#' @param old Logical. Should the old way of detecting outliers be used. Default to FALSE.
#' outlier detection is NOT performed. Useful to get an idea about the initial state of th data.
#' @param normalizeby Should the 2WR matrix be normalized prior to outlier detection, and how.
#' @return A list of class 'phylter' with the 'Initial' (before filtering) and 'Final' (after filtering) states, 
#' or a list of class 'phylterinitial' only, if InitialOnly=TRUE. 
#' @importFrom utils tail combn
#' @importFrom stats hclust as.dist median
#' @importFrom graphics plot
#' @export
phylter<-function(X, bvalue=0, distance="patristic", k=3, k2=k, Norm=TRUE, gene.names=NULL, test.island=TRUE, verbose=TRUE, stop.criteria=1e-5, InitialOnly=FALSE, old=FALSE, normalizeby="row") {
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
	if (is.null(names(X))) X<-rename.genes(X, gene.names=gene.names)
	if (class(X[[1]])=="phylo") matrices <- trees2matrices(X, distance = distance, bvalue = bvalue)
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
	if (Norm==TRUE) matrices<-lapply(matrices, function(x) x/median(x))

#	RES<-DistatisFast(matrices, Norm)
	RES<-DistatisFast(matrices)
	WR<-Dist2WR(RES)
	##qfdkjùjfq
	Initial<-NULL
	Initial$mat.data<-Xsave
	Initial$WR<-WR
	Initial$RV<-RES$RVmat
	Initial$weights<-RES$alpha
	Initial$compromise<-RES$compromise
	Initial$F<-RES$F
	Initial$PartialF<-RES$PartialF
	##New
	Initial$matrices<-matrices

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
			CellsToRemove<-detect.outliers(WR.reordered, k=k,test.island=test.island, old=old, normalizeby=normalizeby)
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
			RES.new<-DistatisFast(matrices.new, Norm)
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
	Final$CELLSREMOVED<-CELLSREMOVED
	Final$Outliers<-CompareBeforeAfter(Xsave, CELLSREMOVED, Final$species.order, which="all")
	Final$CompleteOutliers<-CompareBeforeAfter(Xsave, Final$Outliers, Final$species.order, which="complete")
	#store the way the function was called
	##New 
	Final$matrices<-matrices

	call<-list(call=match.call(), bvalue=bvalue, distance=distance, k=k, Norm=Norm, gene.names=gene.names, test.island=test.island, verbose=verbose, stop.criteria=stop.criteria)

	Result<-list(Initial=Initial, Final=Final, call=call)
	class(Result)<-c("phylter", "list")
	class(Result$Initial)<-c("phylterinitial", "list")
	class(Result$Final)<-c("phylterfinal", "list")

	if (verbose) {
		cat("\n--------\n")
		print(summary(Result))
	}
	return(Result)
}
