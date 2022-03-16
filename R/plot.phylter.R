# Plotting functions for phylter objects. 

#' Plot phylter objects
#' 
#' These functions take objects of class phylter as input
#' and display various plots summarizing the results obtained (see details)
#' 
#' 
#' \itemize{
#'  \item plot(x) and plot.phylter(x) plot the genes found in each species or species 
#' found in each gene as barplots, highlighting the outliers detected.
#'  \item plot2WR(x) plots side by side the initial and the final gene x species (unreable for large datasets) 
#' matrices (the 2WR matrices), highlighting missing data and detected outliers.
#'  \item plotDispersion(x) plots dispersion of data before and after phylter, on a 2D
#' space. Each dot represents a gene-species association. 
#' \item plotRV(x) plots the RV coefficient matrix that descibes all agains all correlations between gene matrices
#' \item plotopti() plots the compromise matrix score at each step of the 
#' optimization.  
#'}
#' @param x The object returned by the 'phylter()' function.
#' @param what Specifies what to plot. If "species", a barplot will show how many
#' genes each species is in, and what proportion of thoses were detected as outliers. 
#' If "genes", a barplot shows how many species contains each gene and how many of 
#' them has been detected as outliers. If "all" (the defaut), the two plots
#' described above are displayed one after the other, prompting the user to press 
#' ENTER to display the next one.  
#' @param layout What layout to use. Do not change if you don't know what it is. 
#' @param sorted Logical. Should bars be sorted by the number of outliers detected. Default
#' to TRUE
#' @param show.missing Logical. Should missing data be represented on the heatmap. If TRUE (the default), white dots show were these missing entries are in both the initial and final 2WR matrices.  
#' @param show.outliers Logical. Should outliers be represented on the heatmap. If TRUE (the default), yellow dots indicate outliers on the final 2WR matrix.
#' @param transpose Logical. If TRUE, the two matrices are piled up instaed of being displayed side by side. Default to FALSE.
#' @param labelnames Logical. If TRUE, the names of labels are indicated on the heatmap. If FALSE they are removed. 
#' This is conveninent when the names of the genes are very long for instance. 
#' @param clust Logical. Should the rows or/and columns of the matrices that are plotted be reorderd
#' prior to plotting. Reordering is based on a hierarchical clustering. Default to FALSE. 
#' This is conveninent when the names of the genes are very long for instance. 

#' @param ... Additional arguments to be passed to plot and print functions.
#' @return The desired plots are returned. Note that you might want to call the pdf(),
#'  png(), jpeg(), or tiff() function first if you want to save the plot(s) to an
#' external file.
#' @examples
#' data(carnivora)
#' 
#' # perform phylter analysis
#' res<-phylter(carnivora)
#' 
#' # plot for each gene the number of outlier species
#' plot(res, "genes")
#' 
#' # plot for each species the number genes where it is outlier
#' plot(res, "species")
#' 
#' # plot the dispersion of data before and after the use of phylter
#' plotDispersion(res)
#' 
#' @importFrom grDevices dev.cur devAskNewPage
#' @importFrom stats relevel
#' @importFrom ggplot2 ggplot aes geom_bar theme element_text labs coord_flip 
#' geom_line geom_point geom_tile scale_color_manual scale_fill_gradient scale_fill_gradient2 element_blank 
#' @importFrom reshape2 melt 
#' @export

plot.phylter<-function(x, what="all", layout=1, sorted=TRUE, ...) {
	## for passing check filters
	namespecies<-NULL
	namegene<-NULL
	number<-NULL
	Species<-NULL
	Genes<-NULL
	##
	sum<-summary(x)
	nbg<-length(sum$initial.nb.sp.per.mat)

	DF_genes<-data.frame(namegene=rep(names(sum$initial.nb.sp.per.mat),2), number=c(sum$nb.sp.removed.per.gene, sum$initial.nb.sp.per.mat-sum$nb.sp.removed.per.gene), Species=factor(c(rep("Removed",nbg), rep("Kept", nbg))))
	DF_genes$Species<-relevel(DF_genes$Species, "Removed")

	if (sorted) {
		GoodOrderGenes<-names(sort(sum$nb.sp.removed.per.gene, decreasing=TRUE))	
		DF_genes$namegene<-factor(DF_genes$namegene, levels=GoodOrderGenes)
	}
	p_genes <- ggplot(DF_genes, aes(x=namegene, y=number, fill=Species)) + geom_bar(stat="identity") + theme(axis.text.x=element_text(angle = 90, hjust = 1, size=2)) + labs(title="Species per gene", x ="Genes", y = "Number of Species") + theme(aspect.ratio=3/6)
	speciesintables<-table(unlist(lapply(x$Initial$mat.data, rownames)))
	nbs<-length(speciesintables)
	speciesinoutliers<-table(x$Final$Outliers[,2])
	speciesinoutliers_reordered<-speciesinoutliers[match(names(speciesintables), names(speciesinoutliers))]
	speciesinoutliers_reordered[is.na(speciesinoutliers_reordered)]<-0
	names(speciesinoutliers_reordered)<-names(speciesintables)
	DF_species<-data.frame(namespecies=rep(names(speciesintables),2), number=c(speciesinoutliers_reordered, speciesintables-speciesinoutliers_reordered), Genes=factor(c(rep("Removed",nbs), rep("Kept", nbs))))
	DF_species$Genes<-relevel(DF_species$Genes, "Removed")

	if (sorted) {
		GoodOrderSpecies<-names(sort(speciesinoutliers_reordered, decreasing=TRUE))
		DF_species$namespecies<-factor(DF_species$namespecies, levels=GoodOrderSpecies)
	}

	p_species <- ggplot(DF_species, aes(x=namespecies, y=number, fill=Genes)) + geom_bar(stat="identity") + theme(axis.text.x=element_text(angle = 90, hjust = 1, size=2), aspect.ratio=3/6) + labs(title="Genes per species",x ="Species", y = "Number of Genes")
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

#' plot2WR
#' @rdname plot.phylter
#' @export

plot2WR<-function(x, show.missing=TRUE, show.outliers=TRUE, transpose=FALSE, clust=FALSE) {
		## for passing check filters
		value<-NULL
		Species<-NULL
		Genes<-NULL
		state<-NULL
		##
		WRi<-x$Initial$WR
		WRf<-x$Final$WR
		sporder<-hclust(dist(WRi))$order
		gnorder<-hclust(dist(t(WRi)))$order
		#we reorder both matrices according to this. 
		if (clust) {
			WRi2<-WRi[sporder, gnorder]
			WRf2<-WRf[sporder, gnorder]
		}
		else {
			WRi2<-WRi
			WRf2<-WRf
		}
		#get xy coordinates of missing data to add to the heatmap
		ExistOrNot<-t(do.call(rbind, lapply(x$Initial$mat.data, function(x,y) match(y, colnames(x)),y=rownames(WRi2))))*0
		ExistOrNot<-ExistOrNot[,gnorder]
		indexOfMissing_i<-which(is.na(ExistOrNot), arr.ind=TRUE)
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
		matoutlier[cbind(match(x$Final$Outliers[,2], rownames(matoutlier)), match(x$Final$Outliers[,1], colnames(matoutlier)))]<-NA
		outliers<-which(is.na(matoutlier), arr.ind=TRUE)
		if (!transpose) outliers[,2]<-outliers[,2]+ncol(WRf2)+1
		outliers<-as.data.frame(outliers)
		outliers$state<-"Discarded outliers"
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
		if (!show.missing & show.outliers) p <- p + geom_point(data=subset(points, state=="Discarded outliers"), aes(x=col, y=row, value=NULL, color=state))
		p <- p+ scale_color_manual(values = c("missing" = "white", "Discarded outliers"="#ffcc00"))
		if (transpose) p <- p + coord_flip()
		p<-p+theme(axis.text=element_text(size=2), aspect.ratio=3/6)
		if (transpose) p <- p + labs(title="gene x species matrix before (bottom) and after (top)\nremoval of phylter-identified outliers")
		else p <- p + labs(title="gene x species matrix before (left) and after (right)\nremoval of phylter-identified outliers")
		print(p)
}

#' plotDispersion
#' @rdname plot.phylter
#' @export


plotDispersion<-function(x) {
	## for passing check filters
	x0<-NULL
	y0<-NULL
	state<-NULL
	##
	step1i<-do.call(rbind, x$Initial$PartialF)
	step2i<-cbind(step1i, x$Initial$F[match(rownames(step1i), rownames( x$Initial$F)),])
	step1f<-do.call(rbind, x$Final$PartialF)
	step2f<-cbind(step1f, x$Final$F[match(rownames(step1f), rownames( x$Final$F)),])
	COO<-rbind(step2i,step2f)
	colnames(COO)<-c("x0","y0","x1","y1")
	COO<-as.data.frame(COO)
	COO$state<-c(rep("Initial", nrow(step2i)), rep("Final", nrow(step2i)))
	COO$genes<-rep(rep(names(x$Final$PartialF), each=length(x$Final$species.order)),2)
	COO$species<-rep(rep(x$Final$species.order,length(x$Final$PartialF)),2)
	COO$outlier<-is.element(paste(COO$species, COO$genes, sep="&&&"), paste(x$Final$Outliers[,2], x$Final$Outliers[,1], sep="&&&")) * (COO$state=="Initial") * 0.25
	rownames(COO)<-NULL
	#center:
	COO$x0<-COO$x0-COO$x1
	COO$x1<-0
	COO$y0<-COO$y0-COO$y1
	COO$y1<-0
	#ggplot(COO, aes(x=x0,y=y0,xend=x1,yend=y1, colour=state)) + geom_segment()
	p <- ggplot(COO, aes(x=x0,y=y0, colour=state)) + geom_point()
	p <- p + labs(caption="One dot = One gene x species association")
	p <- p + labs(title="Dispersion of data before and after removal of phylter-identified outliers")
	print(p)
}


#' plotRV
#' @rdname plot.phylter
#' @export

plotRV<-function(x, what="Initial", labelnames=TRUE, clust=FALSE) {
	## for passing check filters
	Var1<-NULL
	Var2<-NULL
	value<-NULL
	##
	if (what=="Initial") RV<-x$Initial$RV
	if (what=="Final") RV<-x$Final$RV


	if (clust) RV<-RV[hclust(dist(RV))$order,hclust(dist(RV))$order]
	p <- ggplot(melt(RV),aes(x=Var1, y=Var2,fill=value)) + geom_tile() + scale_fill_gradient2(name="RV coefficient", limits=c(-1, 1.01))
	if (!labelnames) {
		p <- p + theme(axis.text.x=element_blank(), axis.text.y=element_blank(), axis.ticks.x=element_blank(), axis.ticks.y=element_blank())
	}
	else {
		p <- p + theme(axis.text.x=element_text(angle = 90, hjust = 1))
	}	

	p <- p + labs(x="Genes",y="Genes", title=paste(what, " vector correlation coefficients (RV) between genes")) + theme(axis.text=element_text(size=2), aspect.ratio=1/1)
	print(p)

#	heatmap(OK$Initial$RV, scale="none", col=magma(100), breaks=seq(-1,1,length.out=101), labCol="", labRow="")
#	require(gplots)
#	heatmap.2(OK$Initial$RV, scale="none", col=magma(100), breaks=seq(-1,1,length.out=101), labCol="", labRow="", density.info="none", trace="none")
}

#' plotopti
#' @rdname plot.phylter
#' @export

plotopti<-function(x) {
	## for passing check filters
	step<-NULL
	score<-NULL
	##

	df<-data.frame(step=1:length(x$Final$AllOptiScores), score=x$Final$AllOptiScores)
	p <- ggplot(df, aes(x=step, y=score)) + geom_line() + geom_point() + labs(title="Evolution of the compromise score with the optimization steps")
	print(p)
}