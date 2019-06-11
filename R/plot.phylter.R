# Plotting functions for phylter objects. 

#' plot.phylter
#' 
#' Plot the initial and final 2WR matrices side by side. Highlights missing data and detected outliers  
#' 
#' @param X The object returned by the 'phylter' function.
#' @param show.missing Logical. Should missing data be represented by white dots? 
#' @param show.outliers Logical. Should outliers be represented by black dots on the final matrix?
#' @param transpose Logical. If FALSE (the default), species are in rows and genes in columns. If TRUE, 
#' rows and columns are inverted.
#' @return A plot of the 2WR matrix.
#' @importFrom grDevices dev.cur devAskNewPage
#' @importFrom utils write.table
#' @importFrom stats relevel
#' @export


plot.phylter<-function(X, what="all", layout=1) {
	sum<-summary(X)
	nbg<-length(sum$initial.nb.sp.per.mat)
	DF_genes<-data.frame(namegene=rep(names(sum$initial.nb.sp.per.mat),2), number=c(sum$nb.sp.removed.per.gene, sum$initial.nb.sp.per.mat-sum$nb.sp.removed.per.gene), Species=c(rep("Removed",nbg), rep("Kept", nbg)))
	DF_genes$Species<-relevel(DF_genes$Species, "Removed")
	p_genes <- ggplot(DF_genes, aes(x=namegene, y=number, fill=Species)) + geom_bar(stat="identity") + theme(axis.text.x=element_text(angle = 90, hjust = 1)) + labs(title="Species per gene", x ="Genes", y = "Number of Species")
	speciesintables<-table(unlist(lapply(X$Initial$mat.data, rownames)))
	nbs<-length(speciesintables)
	speciesinoutliers<-table(X$Final$Outliers[,2])
	speciesinoutliers_reordered<-speciesinoutliers[match(names(speciesintables), names(speciesinoutliers))]
	speciesinoutliers_reordered[is.na(speciesinoutliers_reordered)]<-0
	names(speciesinoutliers_reordered)<-names(speciesintables)
	DF_species<-data.frame(namespecies=rep(names(speciesintables),2), number=c(speciesinoutliers_reordered, speciesintables-speciesinoutliers_reordered), Genes=c(rep("Removed",nbs), rep("Kept", nbs)))
	DF_species$Genes<-relevel(DF_species$Genes, "Removed")
	p_species <- ggplot(DF_species, aes(x=namespecies, y=number, fill=Genes)) + geom_bar(stat="identity") + theme(axis.text.x=element_text(angle = 90, hjust = 1)) + labs(title="Genes per species",x ="Species", y = "Number of Genes")
	if (what=="genes") print(p_genes)
	if (what=="species") print(p_species)
	if (what=="all") {
	    layout(matrix(1:layout, ceiling(sqrt(layout)), byrow = TRUE))
	    if (!devAskNewPage() && !names(dev.cur()) %in% c("pdf", "postscript")) {
	        devAskNewPage(TRUE)
	        on.exit(devAskNewPage(FALSE))
	    }
		print(p_genes)
		print(p_species)
	}

}

plot2WR<-function(X, show.missing=TRUE, show.outliers=TRUE, transpose=FALSE) {
		WRi<-X$Initial$WR
		WRf<-X$Final$WR
		sporder<-hclust(dist(WRi))$order
		gnorder<-hclust(dist(t(WRi)))$order
		#we reorder both matrices according to this. 
		WRi2<-WRi[sporder, gnorder]
		WRf2<-WRf[sporder, gnorder]
		# for (i in 1:nrow(X$Final$Outliers)) {
		# 	WRf2[X$Final$Outliers[i,2],X$Final$Outliers[i,1]]<--diff(range(WRf2, na.rm=TRUE))/100
		# }
		#get xy coordinates of missing data to add to the heatmap
		ExistOrNot<-t(do.call(rbind, lapply(X$Initial$mat.data, function(x,y) match(y, colnames(x)),y=rownames(WRi2))))*0
		ExistOrNot<-ExistOrNot[,gnorder]
		indexOfMissing_i<-which(is.na(ExistOrNot), arr.in=TRUE)
		indexOfMissing_f<-indexOfMissing_i
		if (!transpose) {
			indexOfMissing_f[,2]<-indexOfMissing_f[,2]+ncol(WRf2)+1
		}
		else {
			indexOfMissing_i[,2]<-indexOfMissing_i[,2]+ncol(WRf2)+1
		}		
		if (nrow(indexOfMissing_i)>0) {
			missing<-as.data.frame(rbind(indexOfMissing_i, indexOfMissing_f))
			missing$state<-"missing"
		}
		else {
			missing<-NULL
			show.missing<-FALSE
		}
		#get xy coordinates of outliers
		matoutlier<-WRf2*0
		matoutlier[cbind(match(X$Final$Outliers[,2], rownames(matoutlier)), match(X$Final$Outliers[,1], colnames(matoutlier)))]<-NA
		outliers<-which(is.na(matoutlier), arr.in=TRUE)
		if (!transpose) outliers[,2]<-outliers[,2]+ncol(WRf2)+1
		outliers<-as.data.frame(outliers)
		outliers$state<-"outliers"
		#POINTS
		points<-rbind(missing, outliers)
		#rename columns
		colnames(WRi2)<-paste("Initial - ",colnames(WRi2))
		colnames(WRf2)<-paste("Final - ",colnames(WRf2))
		#make one big matrix:		
		if (transpose) {
			WR<-cbind(WRf2,rep(NA, nrow(WRi2)), WRi2)
		}
		else WR<-cbind(WRi2,rep(NA, nrow(WRi2)), WRf2)

		WRdf<-melt(t(WR))
		colnames(WRdf)<-c("Genes","Species","value")
		p<-ggplot(WRdf, aes(x=Genes, y=Species))
		p <- p + geom_tile(aes(fill = value)) + theme(axis.text.x=element_text(angle = 90, hjust = 1)) 
		#p <- p + scale_fill_gradient(low="#1e6163",high="#ed5f5e", na.value="white") 		
#		p <- p + scale_fill_viridis()
		p <- p + scale_fill_gradient(na.value="white") 		
		if (show.missing & show.outliers) p <- p + geom_point(data=points, aes(x=col, y=row, value=NULL, color=state))
		if (show.missing & !show.outliers) p <- p + geom_point(data=subset(points, state=="missing"), aes(x=col, y=row, value=NULL, color=state))
		if (!show.missing & show.outliers) p <- p + geom_point(data=subset(points, state=="outliers"), aes(x=col, y=row, value=NULL, color=state))
		p <- p+ scale_color_manual(values = c("missing" = "white", "outliers"="#ffcc00"))
		if (transpose) p <- p + coord_flip()
		print(p)
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
	p <- ggplot(COO, aes(x=x0,y=y0, colour=state)) + geom_point()
	p <- p + labs(caption="One gene x species association")
	print(p)
}

##function to write to an easily parsable file the parameters used and the results obtained
writeOutput<-function(X, file="phylter.out") {
	##head
	cat(paste("# \n",sep=""),file=file)
	cat(paste("# -- Phylter v. 0.9 -- \n",sep=""),file=file, append=TRUE)
	cat(paste("# ",date(),"\n", sep=""),file=file, append=TRUE)
	cat(paste("# \n# \n# \n",sep=""),file=file, append=TRUE)
	parameters<-X$call[(names(X$call)!="gene.names")&(names(X$call)!="call")]
	cat(paste("# PARAMETERS\n# \n",sep=""),file=file, append=TRUE)
	cat(paste("# ",paste(names(parameters),parameters, sep="="),sep=""), sep="\n",file=file, append=TRUE)
	phylter_summary<-summary(X)
	cat(paste("# \n# SUMMARY\n# \n",sep=""),file=file, append=TRUE)
	cat(paste("# Total number of outliers detected: ",phylter_summary$nb.outlier.cells,"\n", sep=""),file=file, append=TRUE)
	cat(paste("# Number of complete gene outliers : ",length(phylter_summary$ComplOutGN),"\n", sep=""),file=file, append=TRUE)
	cat(paste("# Number of complete species outliers : ",length(phylter_summary$ComplOutSP),"\n", sep=""),file=file, append=TRUE)
	cat(paste("# Initial score of the compromise: ",X$Final$AllOptiScores[1],"\n", sep=""),file=file, append=TRUE)	
	cat(paste("# Final score of the compromise: ",rev(X$Final$AllOptiScores)[1],"\n", sep=""),file=file, append=TRUE)	
	cat(paste("# Gain: ",phylter_summary$percent.score.increase,"% \n", sep=""),file=file, append=TRUE)
	cat(paste("# Loss (data filtering): ",phylter_summary$percent.data.filtered,"% \n", sep=""),file=file, append=TRUE)
	#	
	if (length(phylter_summary$ComplOutGN)>0) cat(paste("# Outlier gene(s) detected: ",paste(phylter_summary$ComplOutGN, collapse=";"),"\n", sep=""),file=file, append=TRUE)
	if (length(phylter_summary$ComplOutSP)>0) cat(paste("# Outlier species detected: ",paste(phylter_summary$ComplOutSP, collapse=";"),"\n", sep=""),file=file, append=TRUE)
	
	cat(paste("# \n# OUTLIERS (contains complete outliers)\n# \n",sep=""),file=file, append=TRUE)
	cat(paste("# Genes\tSpecies\n",sep=""),file=file, append=TRUE)
	write.table(OK$Final$Outliers, quote=FALSE, col.names=FALSE, row.names=FALSE, sep="\t", file=file, append=TRUE)	
}	

plotRV<-function(X, what="Initial") {
	# if (what=="Initial") {
	# 	ord<-hclust(dist(X$Initial$RV))$order
	# 	RV<-X$Initial$RV[ord,ord]
	# }
	# if (what=="Final") {
	# 	ord<-hclust(dist(X$Final$RV))$order
	# 	RV<-X$Final$RV[ord,ord]
	# }	
	if (what=="Initial") RV<-X$Initial$RV
	if (what=="Final") RV<-X$Final$RV

	RV<-RV[hclust(dist(RV))$order,hclust(dist(RV))$order]
	p <- ggplot(melt(RV),aes(x=Var1, y=Var2,fill=value)) + geom_tile() + scale_fill_gradient2(name="RV coefficient", limits=c(-1, 1.01))
	p <- p + theme(axis.text.x=element_text(angle = 90, hjust = 1)) + labs(x="Genes",y="Genes", title=what)
	print(p)
}

plotopti<-function(X) {
	df<-data.frame(step=1:length(X$Final$AllOptiScores), score=X$Final$AllOptiScores)
	p <- ggplot(df, aes(x=step, y=score)) + geom_line() + geom_point() + labs(title="Evolution of the compromise score with the optimization steps")
	print(p)
}