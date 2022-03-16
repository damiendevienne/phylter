
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
#' write.phylter(res) #writes to the standard output
#' 
#' 
#' @importFrom utils write.table packageVersion
#' @importFrom grDevices dev.off pdf
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
	cat(paste("# Number of genes analyzed: ", dim(ok$Initial$WR)[2],"\n",sep=""),file=file, append=TRUE)
	cat(paste("# Number of species analyzed: ", dim(ok$Initial$WR)[1],"\n",sep=""),file=file, append=TRUE)
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
		pdf(pdfreport.file)
		plot(x, "genes")
		plot(x, "species")
		plotDispersion(x)
		dev.off() #close connection
	}
}	


