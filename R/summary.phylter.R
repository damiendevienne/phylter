summary.phylter<-function(X) {
	percent.score.increase<-round((rev(X$Final$AllOptiScores)[1]-X$Final$AllOptiScores[1])*100, 2)
	nb.outlier.cells<-nrow(X$Final$Outliers)
	initial.nb.sp.per.mat<-unlist(lapply(X$Initial$mat.data, nrow))
	percent.data.filtered<-round(nb.outlier.cells/sum(initial.nb.sp.per.mat)*100,2)
	nb.sp.removed.per.gene<-table(X$Final$Outliers[,1])[names(X$Initial$mat.data)]
	nb.sp.removed.per.gene[is.na(nb.sp.removed.per.gene)]<-0
	names(nb.sp.removed.per.gene)<-names(X$Initial$mat.data)
	print(paste("Gain (concordance between matrices): ",percent.score.increase,"%", sep=""))
	print(paste("Loss (data filtering): ",percent.data.filtered,"%", sep=""))
	#complete outliers: 
	ComplOutSP<-names(which(table(X$Final$Outliers[,2])==length(X$Final$weights)))
	ComplOutGN<-names(which(table(X$Final$Outliers[,1])==length(X$Final$species.order)))

}

print.phylter<-function(X) {
	cat('\n Output of the phylter function. \n\n')
	cat('  $Initial: Initial values before removing outliers:\n')
	cat('    | $mat.data: Initial matrices\n')
	cat('    | $WR: Initial WR matrix\n\n')
	cat('    Final: Final values after removing outliers:\n')
	cat('\n')
}
