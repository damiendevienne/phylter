
#' Write summary of phyter analysis to file(s)
#' 
#' This function a summary of the phylter analysis in a txt file, and can
#' also output a pdf with a full report of the analysis (with graphs).
#' 
#' 
#' @param x The object returned by the 'phylter()' function.
#' @param file Name of the file where to write the text summary of the phylter output.
#' If \code{""} (the default), \code{write.phylyer()} prints to the standard output.
#' @param include.discarded Logical. If TRUE (the default) the elements discarded before the analysis
#' are still in the list of Outliers in the output. Useful for cleaning datasets after a phylter analysis.
#' If \code{""} (the default), \code{write.phylyer()} prints to the standard output.
#' @param pdfreport Logical. Should a full report of the phylter analysis
#  (with graphs) be reported as a pdf file? Default to TRUE
#' @param pdfreport.file If \code{report=TRUE}, name of the pdf file where the 
#' report is written. Default to \code{report.pdf}
#' 
#' @examples
#' data(carnivora) 
#' res<-phylter(carnivora)
#' # write a full report to the standard output
#' write.phylter(res) 
#' 
#' # write a full report to the the file out.txt
#' # write.phylter(res, file="out.txt")
#'  
#' # write a pdf report with all available graphical outputs
#' # write.phylter(res, pdfreport=TRUE) 
#' 
#' @importFrom utils write.table packageVersion
#' @importFrom grDevices dev.off pdf
#' @importFrom graphics plot.new text
#' @export

write.phylter<-function(x, file="", include.discarded=TRUE, pdfreport=FALSE, pdfreport.file="report.pdf") {
	##head
	if (file == "")  {
		file <- stdout()
	}
	cat(paste("# \n",sep=""),file=file)
	cat(paste("# -- Phylter v. ",packageVersion("phylter")," -- \n",sep=""),file=file, append=TRUE)
	cat(paste("# ",date(),"\n", sep=""),file=file, append=TRUE)
	cat(paste("# \n# \n# \n",sep=""),file=file, append=TRUE)
	parameters<-x$call[(names(x$call)!="gene.names")&(names(x$call)!="call")]
	cat(paste("# PARAMETERS\n# \n",sep=""),file=file, append=TRUE)
	cat(paste("# ",paste(names(parameters),parameters, sep="="),sep=""), sep="\n",file=file, append=TRUE)
	phylter_summary<-summary(x)
	cat(paste("# \n# CLEANING STEP\n# \n",sep=""),file=file, append=TRUE)
	cat(paste("# Initial number of genes passed to Phylter: ",phylter_summary$nb.discarded+length(x$Initial$matrices),"\n",sep=""),file=file, append=TRUE)
	cat(paste("# Number of genes discarded before the analysis: ",length(x$DiscardedGenes),"\n",sep=""),file=file, append=TRUE)
	if (length(phylter_summary$nb.discarded)>0) cat(paste("# Genes discarded: ",paste(x$DiscardedGenes, collapse=";"),"\n",sep=""),file=file, append=TRUE)

	cat(paste("# \n# SUMMARY\n# \n",sep=""),file=file, append=TRUE)
	cat(paste("# Number of genes analyzed: ", dim(x$Initial$WR)[2],"\n",sep=""),file=file, append=TRUE)
	cat(paste("# Number of species analyzed: ", dim(x$Initial$WR)[1],"\n",sep=""),file=file, append=TRUE)
	cat(paste("# Total number of outliers detected: ",phylter_summary$nb.outlier.cells,"\n", sep=""),file=file, append=TRUE)
	cat(paste("# Number of complete gene outliers : ",length(phylter_summary$ComplOutGN),"\n", sep=""),file=file, append=TRUE)
	cat(paste("# Number of complete species outliers : ",length(phylter_summary$ComplOutSP),"\n", sep=""),file=file, append=TRUE)
	cat(paste("# Initial score of the compromise: ",x$Final$AllOptiScores[1],"\n", sep=""),file=file, append=TRUE)	
	cat(paste("# Final score of the compromise: ",rev(x$Final$AllOptiScores)[1],"\n", sep=""),file=file, append=TRUE)	
	cat(paste("# Gain: ",phylter_summary$percent.score.increase,"% \n", sep=""),file=file, append=TRUE)
	cat(paste("# Loss (data filtering): ",phylter_summary$percent.data.filtered,"% \n", sep=""),file=file, append=TRUE)
	#	
	if (length(phylter_summary$ComplOutGN)>0) cat(paste("# Outlier gene(s) detected: ",paste(phylter_summary$ComplOutGN, collapse=";"),"\n", sep=""),file=file, append=TRUE)
	if (length(phylter_summary$ComplOutSP)>0) cat(paste("# Outlier species detected: ",paste(phylter_summary$ComplOutSP, collapse=";"),"\n", sep=""),file=file, append=TRUE)
	
	withorwithout<-if (include.discarded) "including" else "excluding"
	cat(paste("# \n# OUTLIERS (",withorwithout," discarded elements, if any)\n# \n",sep=""),file=file, append=TRUE)
	cat(paste("# Genes\tSpecies\n",sep=""),file=file, append=TRUE)
	Out2Write<-if (include.discarded) rbind(x$Final$Outliers, x$Final$Discarded) else x$Final$Outliers
	write.table(Out2Write, quote=FALSE, col.names=FALSE, row.names=FALSE, sep="\t", file=file, append=TRUE)	
	if (pdfreport) {
		f<-pdfreport.file
		if (f=="") f<-"report.pdf"
		#open pdf file
		pdf(pdfreport.file, paper="a4")
		plot.new()
		text1<-paste("- Phylter v. ",packageVersion("phylter")," - ", sep="")
		text1<-paste(text1, date(),sep="\n")
		callstr<-paste(strsplit(toString(x$call),")")[[1]][1],")",sep="")
		text2<-paste("\nCall: ",callstr, "\n\n", sep="")
		text2<-paste(text2, "Initial number of genes passed to Phylter: ",phylter_summary$nb.discarded+length(x$Initial$matrices),"\n",sep="")
		text2<-paste(text2, "Number of genes discarded before the analysis: ",length(x$DiscardedGenes),"\n\n",sep="")
		text2<-paste(text2, "Number of genes analyzed: ", dim(x$Initial$WR)[2],"\n",sep="")
		text2<-paste(text2, "Number of species analyzed: ", dim(x$Initial$WR)[1],"\n",sep="")
		text2<-paste(text2, "Total number of outliers detected: ",phylter_summary$nb.outlier.cells,"\n", sep="")
		text2<-paste(text2, "Number of complete gene outliers : ",length(phylter_summary$ComplOutGN),"\n", sep="")
		text2<-paste(text2, "Number of complete species outliers : ",length(phylter_summary$ComplOutSP),"\n", sep="")
		text2<-paste(text2, "Initial score of the compromise: ",x$Final$AllOptiScores[1],"\n", sep="")	
		text2<-paste(text2, "Final score of the compromise: ",rev(x$Final$AllOptiScores)[1],"\n", sep="")	
		text2<-paste(text2, "Gain: ",phylter_summary$percent.score.increase,"% \n", sep="")
		text2<-paste(text2, "Loss (data filtering): ",phylter_summary$percent.data.filtered,"% \n", sep="")
		text3<-"See text report for full details"
		text(0.5,0.9,text1)		
		text(0,0.5,text2, cex=0.8, font=2, adj=0, col="#333333")		
		text(0.5,0,text3, cex=0.8)		
		plot(x, "genes")
		plot(x, "species")
		plotDispersion(x)
		plot2WR(x)
		plotRV(x, "Initial")
		plotRV(x, "Final")
		plotopti(x)
		dev.off() #close connection
	}
}	


