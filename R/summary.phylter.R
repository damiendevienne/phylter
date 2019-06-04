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
	ComplOutSP<-names(which(table(X$Final$Outliers[,2])==length(X$Final$weights)))
	ComplOutGN<-names(which(table(X$Final$Outliers[,1])==length(X$Final$species.order)))
	#result
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

print.summary.phylter<-function(x) {
	if (!inherits(x, "summary.phylter"))
		stop("'x' must inherit from class summar.phylter")
	cat(paste("Number of outlier cells detected: ",x$nb.outlier.cells,"\n", sep=""))
	cat(paste("Gain (concordance between matrices): ",x$percent.score.increase,"% \n", sep=""))
	cat(paste("Loss (data filtering): ",x$percent.data.filtered,"% \n", sep=""))
}

print.phylter<-function(X) {
	cat('\n Output of the phylter function. \n\n')
	cat('  $Initial: Initial values before removing outliers:\n')
	cat('    | $mat.data: Initial matrices\n')
	cat('    | $WR: Initial WR matrix\n\n')
	cat('    Final: Final values after removing outliers:\n')
	cat('\n')
}
