rm(list=ls())
tre<-read.tree("datasets/droso.tre")
clean<-function(tr) {
	n<-unlist(lapply(strsplit(tr$tip.label, "_"), function(x) x[[1]]))
	tr$tip.label<-n
	if (sum(duplicated(n))>0) {
		tr$ismono<-is.monophyletic(tr, which(duplicated(n)))
		tr<-drop.tip(tr, which(duplicated(n)))
	}
	else tr$ismono<-TRUE
	tr
}
trees<-lapply(tre, clean)
ismono<-unlist(lapply(trees, function(x) x$ismono)) #if true, in these trees the DROME, if in more than 1 copy, are monophyletic
trees<-trees[ismono]
