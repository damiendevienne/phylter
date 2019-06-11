# Summary and print functions for objects of class phylter, phylterfinal and phylterinitial 
# No need to document those?

#' summary.phylter
#' 
#' summary.phylter TODO.
#' 
#' @param X Object returned by function 'phylter()'.
#' @return Print formatting
#' @export

summary.phylter<-function(X) {
	res<-NULL
	percent.score.increase<-round((rev(X$Final$AllOptiScores)[1]-X$Final$AllOptiScores[1])*100, 2)
	nb.outlier.cells<-nrow(X$Final$Outliers)
	initial.nb.sp.per.mat<-unlist(lapply(X$Initial$mat.data, nrow))
	percent.data.filtered<-round(nb.outlier.cells/sum(initial.nb.sp.per.mat)*100,2)
	nb.sp.removed.per.gene<-table(X$Final$Outliers[,1])[names(X$Initial$mat.data)]
	nb.sp.removed.per.gene[is.na(nb.sp.removed.per.gene)]<-0
	names(nb.sp.removed.per.gene)<-names(X$Initial$mat.data)
	#complete outliers: 
	ComplOutSP<-X$Final$CompleteOutliers$ComplOutSP
	ComplOutGN<-X$Final$CompleteOutliers$ComplOutGN
	#result
	res$nb.outlier.cells<-nb.outlier.cells
	res$percent.score.increase<-percent.score.increase
	res$percent.data.filtered<-percent.data.filtered
	res$initial.nb.sp.per.mat<-initial.nb.sp.per.mat
	res$nb.sp.removed.per.gene<-nb.sp.removed.per.gene
	res$ComplOutGN<-ComplOutGN
	res$ComplOutSP<-ComplOutSP
	res$keep.species<-X$call$keep.species
	class(res)<-"summary.phylter"
	res
}

#' print.summary.phylter
#' 
#' print.summary.phylter TODO.
#' 
#' @param X Object returned by function 'phylter()'.
#' @return Print formatting
#' @export

print.summary.phylter<-function(x) {
	if (!inherits(x, "summary.phylter"))
		stop("'x' must inherit from class summar.phylter")
	cat("\n")
	cat(paste("Total number of outliers detected: ",x$nb.outlier.cells,"\n", sep=""))
	cat(paste("  Number of complete gene outliers : ",length(x$ComplOutGN),"\n", sep=""))
	cat(paste("  Number of complete species outliers : ",ifelse(x$keep.species==TRUE, " NA*\n  *set keep.species=FALSE to allow this detection",length(x$ComplOutSP)),"\n", sep=""))
	cat("\n")
	cat(paste("Gain (concordance between matrices): ",x$percent.score.increase,"% \n", sep=""))
	cat(paste("Loss (data filtering): ",x$percent.data.filtered,"% \n", sep=""))
}

print.phylter<-function(X) {
	cat("Phylter Analysis\nList of class phylter\n\n", sep="")
	cat("Call: ")
	print(X$call$call)
	cat("\n")
	cat('$Initial\tInitial matrices and values, before optimization\n')
	cat('$Final\t\tFinal matrices, scores, outliers, after optimization\n\n')
	cat("\n\nTips:\n")
	cat('   Use summary(X) to get an overview of the results.\n')
	cat('   Use plot(X) to see the distribution of outliers\n')
	cat('   Use plot2WR(X) to compare WR matrices before and after\n')
	cat('   Use write.phylter(X) to write the results to an easily parsable file\n')
	cat('   Use plotDispersion(X) to compare Distatis projections before and after\n')	
}

print.phylterinitial<-function(X) {
	cat("Phylter Analysis - initial state\nList of class phylterinitial\n\n", sep="")

	Object<-paste("$",names(X), sep="")
	Dimension<-unlist(lapply(X, function(x) ifelse(class(x)=="matrix",paste(dim(x),collapse=" x "), paste(length(x)))))
	Content<-c(
		"List of original distance matrices, one per gene",
		"Species x Genes reference matrix",
		"Genes x Genes RV correlation coefficients matrix",
		"Weight of each gene in the compromise",
		"Species x Species compromise matrix",
		"Distatis coordinates of compromise",
		"Distatis coordinates of gene matrices (list)"
		)
	DF<-data.frame(Object, Dimension, Content, row.names=NULL)
	print(DF, right=FALSE)
}

print.phylterfinal<-function(X) {
	cat("Phylter Analysis - final state\nList of class phylterfinal\n\n", sep="")

	Object<-paste("$",names(X), sep="")
	Dimension<-unlist(lapply(X, function(x) ifelse(class(x)=="matrix",paste(dim(x),collapse=" x "), paste(length(x)))))
	Content<-c(
		"Species x Genes reference matrix",
		"Genes x Genes RV correlation coefficients matrix",
		"Weight of each gene in the compromise",
		"Species x Species compromise matrix",
		"Distatis coordinates of compromise",
		"Distatis coordinates of gene matrices (list)",
		"Name and order of species",
		paste("Evolution of quality of compromise (",length(X$AllOptiScores)," steps)",sep=""),
		"Index of cells removed (may contain imputed cells)",
		"Outliers detected (one row = one outlier cell)",
		"Complete outliers (Gene and Species, if any)"
		)
	DF<-data.frame(Object, Dimension, Content, row.names=NULL)
	print(DF, right=FALSE)
}

