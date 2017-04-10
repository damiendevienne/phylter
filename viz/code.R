require(ape)
require(phangorn)
require(ade4)
require(scales)

t<-read.tree("Aguileta-et-al-2008_TREES.txt")
tc<-lapply(t, getClans)
## we remove clades with only one 0 or one 1
tc<-lapply(tc, function(x) x[(apply(x,1,sum)!=1)&(apply(x,1,sum)!=(ncol(x)-1)),])



nam<-colnames(tc[[1]])


TC<-lapply(tc, function(x,y) apply(x[,y],1,paste, collapse=""), y=nam)



first<-unique(unlist(TC))
##and we remove the copies (those that have one instead of zeros)
revert<-gsub("2","0",gsub("0","1",gsub("1","2",first)))
KK<-sapply(first, function(x,y) which(y==x), y=revert)
##we remove those that have a number in KK smaller than their rank.
UNIKBP<-names(KK[KK>1:length(KK)])

TCf<-lapply(TC, function(x, y) x[is.element(x,y)], y=UNIKBP)

###WE NEED to decide how we place the different bipartitions around a circle.

circ<-function(n) {
    alpha0 = pi/2
    alpha <- alpha0 - (1:n) * 2 * pi/n
    x <- cos(alpha)
    y <- sin(alpha)
    return(cbind(x,y))
}

ok<-circ(length(UNIKBP))
rownames(ok)<-sample(names(sort(table(unlist(TCf)))))

plot(ok, pch=19, xlab="", ylab="", frame=FALSE, axes=FALSE)
for (i in 1:length(TCf)) {
    where<-is.element(rownames(ok), TCf[[i]])
    polygon(ok[where,], lwd=2, border=alpha("blue", 0.1))
}









###OTHER approach: one circle per species and per gene
tr<-read.tree("Aguileta-et-al-2008_TREES.txt")
tr<-lapply(tr, compute.brlen, 1)
TAB<-lapply(tr, cophenetic)
nam<-tr[[1]]$tip.label
TAB<-lapply(TAB, function(x,y) x[y,y],y=nam)






##let's select one species (SP)
par(mfrow=c(5,5))
par(mar=c(0,0,0,0))
par(oma=c(0,0,0,0))
for (j in 1:length(nam)) { ##for each species
    SP<-nam[j]
    T1<-lapply(TAB, function(x) (x[SP,nam]))
    T1m<-matrix(unlist(T1), nrow=length(tr), byrow=TRUE)
    ##T1m gives 1 plot corresponding to "Kla" for each gene.
    Means.T1m<-apply(T1m, 2, mean)
    
    ##CIRCLE:
    xc<-rep(1,length(nam)+1)*cos(seq(0,2*pi,length.out=length(nam)+1))
    yc<-rep(1,length(nam)+1)*sin(seq(0,2*pi,length.out=length(nam)+1))
    ##we check angles
    alphas<-seq(0,2*pi,length.out=length(nam)+1)
    alphas<-alphas[1:length(nam)]

    xc<-xc[1:21]
    yc<-yc[1:21]
    alphas<-alphas[1:21]
    ##for each gene, the ray is given by the proportion:
    
    GENEi<-NULL
    plot(4*xc,4*yc,type="n", xlim=c(-4,4), ylim=c(-4,4), frame.plot=FALSE, axes=FALSE, xlab="", ylab="")
    text(4*xc,4*yc,labels=nam, col="light grey")
    for (i in 1:246) {
        genei<-T1m[i,]/Means.T1m
        ##NEW! 
#        genei<-1+abs(1-genei)
        genei[is.na(genei)]<-1
        GENEi<-c(GENEi, genei)
        x<-genei*cos(alphas)
        y<-genei*sin(alphas)
        x[is.na(x)]<-0
        y[is.na(y)]<-0
        ##    plot(x,y,type="n", xlim=c(-2,2), ylim=c(-2,2), frame.plot=FALSE, axes=FALSE, xlab="", ylab="")
        polygon(xc,yc, border="light grey", lwd=0.54)
        polygon(x,y,border="red", lwd=0.8)
        text(-3.5,-3.5,SP,cex=2)
    }
}



#####WE CAN ALSO DO THAT BY GENE

par(mfrow=c(10,10))
par(mar=c(0,0,0,0))
par(oma=c(0,0,0,0))
for (i in 201:246) { ##25 first genes
    plot(x,y,type="n", xlim=c(-4,4), ylim=c(-4,4), frame.plot=FALSE, axes=FALSE, xlab="", ylab="")
    for (j in 1:length(nam)) { ##for each speciew
        SP<-nam[j]
        T1<-lapply(TAB, function(x,y) (x[SP,nam]))
        T1m<-matrix(unlist(T1), nrow=length(tr), byrow=TRUE)
        ##T1m gives 1 plot corresponding to "Kla" for each gene.
        Means.T1m<-apply(T1m, 2, mean)
        genei<-T1m[i,]/Means.T1m
        ##        genei<-1+abs(1-genei)

        x<-genei*cos(alphas)
        y<-genei*sin(alphas)
        x[is.na(x)]<-0
        y[is.na(y)]<-0
        ##    plot(x,y,type="n", xlim=c(-2,2), ylim=c(-2,2), frame.plot=FALSE, axes=FALSE, xlab="", ylab="")
        polygon(xc,yc, border="light grey", lwd=0.5)
        points(x,y,pch=19, cex=0.2, col="red")
        polygon(x,y,border="red", lwd=0.1)
#        text(-2.5,-2.5,paste("gene",i,sep=" "),cex=1)
    }
}









##BY SPECIES AGAIN
for (j in 1:21) {
     SP<-nam[j]
     T1<-lapply(TAB, function(x,y) (x[SP,nam]))
     T1m<-matrix(unlist(T1), nrow=length(tr), byrow=TRUE)
    for (i in 1:246) {
        ii<-apply(T1m, 1, function(x) cor.test(x, T1m[i,])$estimate)
    }
}
