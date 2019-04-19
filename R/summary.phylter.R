summary.phylter<-function(X) {
	percent.score.increase<-round((rev(X$Final$AllOptiScores)[1]-X$Final$AllOptiScores[1])*100, 2)
	nb.outlier.cells<-nrow(X$Final$Outliers)
	initial.nb.sp.per.mat<-unlist(lapply(X$Initial$mat.data, nrow))
	percent.data.filtered<-round(nb.outlier.cells/sum(initial.nb.sp.per.mat)*100,2)
}