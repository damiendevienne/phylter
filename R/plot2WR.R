# plot the 2WR matrices, either the initial one of the final one..

#' plot2WR
#' 
#' plot the 2WR matrices, either the Initial one or the Final one.  
#' 
#' @param X The object returned by the 'phylter' function.
#' @param what Should the "Initial" or "Final" 2WR matrix be plotted?
#' @return A plot of the 2WR matrix.
#' @export

plot2WR<-function(X) {
		WRi<-X$Initial$WR
		WRf<-X$Final$WR
		sporder<-hclust(dist(WRi))$order
		gnorder<-hclust(dist(t(WRi)))$order
		#we reorder both matrices according to this. 
		WRi2<-WRi[sporder, gnorder]
		WRf2<-WRf[sporder, gnorder]
		for (i in 1:nrow(X$Final$Outliers)) {
			WRf2[X$Final$Outliers[i,2],X$Final$Outliers[i,1]]<--diff(range(WRf2))/100
		}
		colnames(WRi2)<-paste("Initial - ",colnames(WRi2))
		colnames(WRf2)<-paste("Final - ",colnames(WRf2))
		#make one big matrix:
		WR<-cbind(WRi2,rep(NA, nrow(WRi2)), WRf2)
		WRdf<-melt(t(WR))
		colnames(WRdf)<-c("Genes","Species","value")

		ggplot(WRdf, aes(x=Genes, y=Species, fill=value)) + geom_tile() + theme(axis.text.x=element_text(angle = 90, hjust = 1)) + scale_fill_gradientn(colors=c("grey",heat.colors(100)), na.value="white") 
}

plotDispersion<-function(X) {
	step1i<-do.call(rbind, X$Initial$PartialF)
	step2i<-cbind(step1i, X$Initial$F[match(rownames(step1i), rownames( X$Initial$F)),])
	step1f<-do.call(rbind, X$Final$PartialF)
	step2f<-cbind(step1f, X$Final$F[match(rownames(step1f), rownames( X$Final$F)),])
	COO<-rbind(step2i,step2f)
	colnames(COO)<-c("x0","y0","x1","y1")
	COO<-as.data.frame(COO)
	COO$state<-c(rep("Initial", nrow(step2i)), rep("Final", nrow(step2i)))
	COO$genes<-rep(rep(names(X$Final$PartialF), each=length(X$Final$species.order)),2)
	COO$species<-rep(rep(X$Final$species.order,length(X$Final$PartialF)),2)
	COO$outlier<-is.element(paste(COO$species, COO$genes, sep="&&&"), paste(X$Final$Outliers[,2], X$Final$Outliers[,1], sep="&&&")) * (COO$state=="Initial") * 0.25
	rownames(COO)<-NULL
	#center:
	COO$x0<-COO$x0-COO$x1
	COO$x1<-0
	COO$y0<-COO$y0-COO$y1
	COO$y1<-0
	#ggplot(COO, aes(x=x0,y=y0,xend=x1,yend=y1, colour=state)) + geom_segment()
	ggplot(COO, aes(x=x0,y=y0, colour=state)) + geom_point()
}

##function to write to an easily parsable file the parameters used and the results obtained
writeOutput<-function(X, file="phylter.out") {
	##head
	cat(paste("# \n",sep=""),file=file)
	cat(paste("# -- Phylter v. 0.9 -- \n",sep=""),file=file, append=TRUE)
	cat(paste("# \n# \n# \n",sep=""),file=file, append=TRUE)
	##
	parameters<-X$call[names(X$call)!="gene.names"]
	cat(paste("# PARAMETERS\n# \n",sep=""),file=file, append=TRUE)
	cat(paste("# ",paste(names(parameters),parameters, sep="="),sep=""), sep="\n",file=file, append=TRUE)
	phylter_summary<-summary(X)
	cat(paste("# \n# SUMMARY\n# \n",sep=""),file=file, append=TRUE)
	cat(paste("# Total number of outliers detected: ",phylter_summary$nb.outlier.cells,"\n", sep=""),file=file, append=TRUE)
	cat(paste("# Number of complete gene outliers : ",length(phylter_summary$ComplOutGN),"\n", sep=""),file=file, append=TRUE)
	cat(paste("# Number of complete species outliers : ",length(phylter_summary$ComplOutSP),"\n", sep=""),file=file, append=TRUE)
	cat("# \n",file=file, append=TRUE)
	cat(paste("# Initial score of the compromise: ",X$Final$AllOptiScores[1],"\n", sep=""),file=file, append=TRUE)	
	cat(paste("# Final score of the compromise: ",rev(X$Final$AllOptiScores)[1],"\n", sep=""),file=file, append=TRUE)	
	cat(paste("# Gain: ",phylter_summary$percent.score.increase,"% \n", sep=""),file=file, append=TRUE)
	cat("# \n",file=file, append=TRUE)
	
	cat(paste("# Loss (data filtering): ",phylter_summary$percent.data.filtered,"% \n", sep=""),file=file, append=TRUE)
	cat(paste("# \n# OUTLIERS\n# \n",sep=""),file=file, append=TRUE)
	cat(paste("# Complete gene outliers:\n# ",phylter_summary$ComplOutGN,"\n", sep=""),file=file, append=TRUE)
	cat(paste("# Complete species outliers:\n# ",phylter_summary$ComplOutSP,"\n", sep=""),file=file, append=TRUE)
	cat(paste("# \n# ALL OUTLIERS\n# \n",sep=""),file=file, append=TRUE)

}	