#' Get summary for phylter objects
#' 
#' Get the summary of an object of class \code{phylter}
#' 
#' @param object Object returned by function 'phylter()'.
#' @param ... Additional arguments affecting the summary produced.
#' @return On object of class \code{summmary.phylter}
#' @examples 
#' data(carnivora)
#' res<-phylter(carnivora)
#' summary<-summary(res)
#' class(summary)
#' @export

summary.phylter<-function(object, ...) {
	res<-NULL
	percent.score.increase<-round((rev(object$Final$AllOptiScores)[1]-object$Final$AllOptiScores[1])*100, 2)
	nb.discarded<-length(object$DiscardedGenes)
	nb.outlier.cells<-nrow(object$Final$Outliers)
	initial.nb.sp.per.mat<-unlist(lapply(object$Initial$mat.data, nrow))
	percent.data.filtered<-round(nb.outlier.cells/sum(initial.nb.sp.per.mat)*100,2)
	nb.sp.removed.per.gene<-table(object$Final$Outliers[,1])[names(object$Initial$mat.data)]
	nb.sp.removed.per.gene[is.na(nb.sp.removed.per.gene)]<-0
	names(nb.sp.removed.per.gene)<-names(object$Initial$mat.data)
	#complete outliers: 
	ComplOutSP<-object$Final$CompleteOutliers$ComplOutSP
	ComplOutGN<-object$Final$CompleteOutliers$ComplOutGN
	#result
	res$nb.discarded<-nb.discarded
	res$nb.outlier.cells<-nb.outlier.cells
	res$percent.score.increase<-percent.score.increase
	res$percent.data.filtered<-percent.data.filtered
	res$initial.nb.sp.per.mat<-initial.nb.sp.per.mat
	res$nb.sp.removed.per.gene<-nb.sp.removed.per.gene
	res$ComplOutGN<-ComplOutGN
	res$ComplOutSP<-ComplOutSP
	class(res)<-"summary.phylter"
	res
}

#' print summary of phylter objects
#' 
#' Prints on screen the summary of an object of class summmary.phylter as returned by the \code{summary.phylter()} function.
#' 
#' @param x Object returned by function 'summary.phylter()'.
#' @param ... Additional arguments.
#' @return NA   
#' @examples
#' data(carnivora)
#' res<-phylter(carnivora)
#' summary<-summary(res)
#' print(summary)
#' @export

print.summary.phylter<-function(x, ...) {
	if (!inherits(x, "summary.phylter"))
		stop("'x' must inherit from class summar.phylter")
	cat("\n")
	cat(paste("Total number of outliers detected: ",x$nb.outlier.cells,"\n", sep=""))
	cat(paste("  Number of complete gene outliers : ",length(x$ComplOutGN),"\n", sep=""))
	cat(paste("  Number of complete species outliers : ",length(x$ComplOutSP),"\n", sep=""))
	cat("\n")
	cat(paste("Gain (concordance between matrices): ",x$percent.score.increase,"% \n", sep=""))
	cat(paste("Loss (data filtering): ",x$percent.data.filtered,"% \n", sep=""))
}


#' Print phylter objects
#' 
#' Print on screen a simple description of the content of phylter objects
#' 
#' @param x Object returned by function 'phylter()'.
#' @param ... Additional arguments.
#' @return Print formatting   
#' @examples
#' data(carnivora)
#' res<-phylter(carnivora)
#' print(res)
#' @export


print.phylter<-function(x, ...) {
	cat("Phylter Analysis\nList of class phylter\n\n", sep="")
	cat("Call: ")
	print(x$call$call)
	cat("\n")
	cat('$Initial\tInitial matrices and values, before optimization\n')
	cat('$Final\t\tFinal matrices, scores, outliers, after optimization\n')
	cat('$DiscardedGenes\tList of discarded genes (not analyzed by phylter)\n\n')	
	cat("\n\nTips:\n")
	cat('   Use summary(x) to get an overview of the results.\n')
	cat('   Use plot(x) to see the distribution of outliers\n')
	cat('   Use plot2WR(x) to compare WR matrices before and after\n')
	cat('   Use write.phylter(x) to write the results to an easily parsable file\n')
	cat('   Use plotDispersion(x) to compare Distatis projections before and after\n')	
}


#' print objects of class phylterinitial
#' 
#' Print on screen a simple description of the content of objects of class phylterinitial
#' 
#' @param x Object present in the \code{$initial} element of the object returned by function \code{phylter()}.
#' @param ... Additional arguments.
#' @return NA 
#' @examples
#' data(carnivora)
#' res<-phylter(carnivora)
#' print(res$Initial)
#' @export


print.phylterinitial<-function(x, ...) {
	cat("Phylter Analysis - initial state\nList of class phylterinitial\n\n", sep="")

	Object<-paste("$",names(x), sep="")
	Dimension<-unlist(lapply(x, function(x) ifelse(class(x)[1]=="matrix",paste(dim(x),collapse=" x "), paste(length(x)))))
	Content<-c(
		"List of original distance matrices, one per gene",
		"Species x Genes reference matrix",
		"Genes x Genes RV correlation coefficients matrix",
		"Weight of each gene in the compromise",
		"Species x Species compromise matrix",
		"Distatis coordinates of compromise",
		"Distatis coordinates of gene matrices (list)", 
		"Species x Species gene matrices (list)"
		)
	DF<-data.frame(Object, Dimension, Content, row.names=NULL)
	print(DF, right=FALSE)
}



#' print objects of class phylterfinal
#' 
#' Print on screen a simple description of the content of objects of class phylterfinal
#' 
#' @param x Object returned by function 'phylter()'.
#' @param ... Additional arguments.
#' @return NA   
#' @examples
#' data(carnivora)
#' res<-phylter(carnivora)
#' print(res$Final)
#' @export


print.phylterfinal<-function(x, ...) {
	cat("Phylter Analysis - final state\nList of class phylterfinal\n\n", sep="")

	Object<-paste("$",names(x), sep="")
	Dimension<-unlist(lapply(x, function(x) ifelse(class(x)[1]=="matrix",paste(dim(x),collapse=" x "), paste(length(x)))))
	Content<-c(
		"Species x Genes reference matrix",
		"Genes x Genes RV correlation coefficients matrix",
		"Weight of each gene in the compromise",
		"Species x Species compromise matrix",
		"Distatis coordinates of compromise",
		"Distatis coordinates of gene matrices (list)",
		"Name and order of species",
		paste("Evolution of quality of compromise (",length(x$AllOptiScores)," steps)",sep=""),
		"Index of cells removed (may contain imputed cells)",
		"Outliers detected (one row = one outlier cell)",
		"Complete outliers (Gene and Species, if any)",
		"Species x Species gene matrices (list)"
		)
	if (!is.null(x$Discarded)) Content<-c(Content, "Discarded elements (one row = one species/gene discarded)")
	DF<-data.frame(Object, Dimension, Content, row.names=NULL)
	print(DF, right=FALSE)
}

