
#' Write phylter summary to file
#' 
#' This function writes to a file a summary of the phylter analysis.
#' 
#' 
#' @param x The object returned by the 'phylter()' function.
#' @param file Name of the file where to write the summary of the phylter output.
#' @importFrom utils write.table
#' @export

write.phylter<-function(x, file="phylter.out") {
	##head
	cat(paste("# \n",sep=""),file=file)
	cat(paste("# -- Phylter v. 0.9 -- \n",sep=""),file=file, append=TRUE)
	cat(paste("# ",date(),"\n", sep=""),file=file, append=TRUE)
	cat(paste("# \n# \n# \n",sep=""),file=file, append=TRUE)
	parameters<-x$call[(names(x$call)!="gene.names")&(names(x$call)!="call")]
	cat(paste("# PARAMETERS\n# \n",sep=""),file=file, append=TRUE)
	cat(paste("# ",paste(names(parameters),parameters, sep="="),sep=""), sep="\n",file=file, append=TRUE)
	phylter_summary<-summary(x)
	cat(paste("# \n# SUMMARY\n# \n",sep=""),file=file, append=TRUE)
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
	
	cat(paste("# \n# OUTLIERS (contains complete outliers)\n# \n",sep=""),file=file, append=TRUE)
	cat(paste("# Genes\tSpecies\n",sep=""),file=file, append=TRUE)
	write.table(x$Final$Outliers, quote=FALSE, col.names=FALSE, row.names=FALSE, sep="\t", file=file, append=TRUE)	
}	