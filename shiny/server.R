
server <-function(input,output,session){

 ##reactives
  outVar <- reactive({
    trees = trees()
    sps = trees[[1]]$tip.label
    if (!is.null(sps)){
      vars <- all.vars(parse(text = sps))
    }
  })
  outliers <- reactive({
    input$Outlier
  })
  RES <- reactive({
    trees = trees()
    if(is.null(names(trees))){
      names(trees)=as.character(c(1:length(trees)))
    }
    Phylter(trees, distance = input$choice, k=input$k, thres = 0.5)
  })

  dist <- reactive({
    trees = trees()
    if(is.null(names(trees))){
      names(trees)=as.character(c(1:length(trees)))
    }
    mat <- trees2matrices.Distatis(trees, input$choice)
    mat2Dist(mat)
  })
  choice <- reactive({
    input$choice
  })
  TreesList <- reactive({
    input$Treeslist
  })
  
  WR <- reactive({
    Dist2WR(dist())
  })
  trees <- reactive({
    trees = read.tree(input$trees$datapath, keep.multi = TRUE)
  })
  ranges <- reactiveValues(x = NULL, y = NULL)
  
###########-----------------------------------
  
  #quand on ajoute un arbre
  observeEvent(input$trees, {
    trees = trees()
    if(is.null(names(trees))){
      names(trees)=as.character(c(1:length(trees)))
    }
    FS=dist()$res4Splus$F
    PF=dist()$res4Splus$PartialF
    C = dist()$res4Cmat$C
    if(is.null(dimnames(PF)[[3]])){
      if(!is.null(names(trees))){
        dimnames(PF)[[3]] = names(trees)
      }
      else {
        dimnames(PF)[[3]] = as.character(c(1:length(trees)))
      }
    }
    G=dist()$res4Cmat$G
    if(is.null(rownames(G))){
      if(!is.null(names(trees))){
        rownames(G)<-names(trees)
      }
      else{
        rownames(G)<-as.character(c(1:length(trees)))
      }
    }
    ##Visualize trees
    observeEvent(input$SumitTrees,{
      List = TreesList()
      S = strsplit(List, ",")
      S = S[[1]]
      output$plot1 <- renderPlot({
        if (length(S)==1){par(mfrow = c(1,1))}
        if (length(S)>=2 && length(S)<=8){par(mfrow = c(ceiling(length(S)/2),2))}
        if (length(S)>8 && length(S)<=30){par(mfrow = c(ceiling(length(S)/5),5))}
        if (length(S)>30){par(mfrow = c(ceiling(length(S)/8),8))}
        par(mar = c(0,0,0,0))
        for (i in 1:length(S)){
          ape::plot.phylo(trees[[S[i]]])
        }
      })
    })
    
    ##Update Genetrees
    n=names(trees)
    c=1
    l1 = 20
    Partition = list()
    while (l1 < length(n)){
      p=n[l1]
      k=paste(c, "-", p, sep="")
      Partition = append(Partition, k)
      l1 =l1+20
      c=c+20
    }
    if (l1 >= length(n)){
      k=paste(c, "-", length(n), sep="")
      Partition = append(Partition, k)
    }
    observe({
      updateSelectInput(session, "Geneslist", choices = Partition)
    })
    
    observeEvent(outliers(), {
      ##event = DETECTION OUTLIER ON
      if (outliers()== "on"){
        output$outputPhylter1 <- renderText ({
          if(length(RES()$Complete$outgn)!=0){
            gn=paste(RES()$Complete$outgn, ";")
            append (gn, "outlier genes: ", after=0)
          }
          else{
            "no outlier gene detected"
          }
        })
        output$outputPhylter2 <- renderText ({  
          if(length(RES()$Complete$outsp)!=0){
            sp=paste(RES()$Complete$outsp, ";")
            append (sp, "outlier species: ", after=0)
          }
        else{
          "no outlier specie detected"
        }
      }) 
      output$outputPhylter3 <- renderText ({  
        if(!is.null(RES()$CellByCell$outcell)){
          ce = paste(RES()$CellByCell$outcell[,1],RES()$CellByCell$outcell[,2], ";")
          append (ce, "outlier gène/species: ", after=0)
        }
        else{
          "no outlier cell detected"
        }
      })
      trees = trees()
      G=dist()$res4Cmat$G
      if(is.null(rownames(G))){
        if(!is.null(names(trees))){
          rownames(G)<-names(trees)
        }
        else{
          rownames(G)<-as.character(c(1:length(trees)))
        }
      }
      output$plot3 <- renderPlot({
        gn=paste(RES()$Complete$outgn)
        participant.colors <- rownames(G)
        if (length(RES()$Complete$outgn)!=0){
          for (j in 1:length(gn)){
            participant.colors[which(participant.colors==gn[j])] = "red1"
          }
          participant.colors[which(participant.colors!="red1")] = "grey50"
        }
        else{
          participant.colors = "grey50"
        }
        plot(G, col = participant.colors, cex = 3.5)
        text(G, labels = rownames(G))
      })
      ##Distance by genes
      output$plot7 <- renderPlot({
        Geneslist = input$Geneslist
        Geneslist = strsplit(Geneslist, "-")
        TAB<-trees2matrices.Distatis(trees, distance=choice())
        nam<-trees[[1]]$tip.label
        #TAB<-lapply(TAB, function(x,y) x[y,y],y=nam)
        par(mfrow=c(5,4))
        par(mar=c(1,1,1,1))
        par(oma=c(0,0,0,0))
        i1 = as.integer(Geneslist[[1]][1])
        i2 = as.integer(Geneslist[[1]][2])
        for (i in i1:i2) {
          gn=paste(RES()$Complete$outgn)
          if (i%in%gn){
            colo = "red"
          }
          else{
            colo = "black"
          }
          plot(0,0,type="n", xlim=c(-4,4), ylim=c(-4,4), frame.plot=FALSE, axes=FALSE, xlab="", ylab="", main=i, col.main=colo, cex.main = 1.5)
          for (j in 1:length(nam)) { ##for each speciew
            SP<-nam[j]
            T1<-lapply(TAB, function(x,y) (x[SP,nam]))
            T1m<-matrix(unlist(T1), nrow=length(trees), byrow=TRUE)
            Means.T1m<-apply(T1m, 2, mean)
            genei<-T1m[i,]/Means.T1m
            xc<-rep(1,length(nam)+1)*cos(seq(0,2*pi,length.out=length(nam)+1))
            yc<-rep(1,length(nam)+1)*sin(seq(0,2*pi,length.out=length(nam)+1))
            ##we check angles
            alphas<-seq(0,2*pi,length.out=length(nam)+1)
            alphas<-alphas[1:length(nam)]
            x<-genei*cos(alphas)
            y<-genei*sin(alphas)
            x[is.na(x)]<-0
            y[is.na(y)]<-0
            polygon(xc,yc, border="light grey", lwd=0.5)
            points(x,y,pch=19, cex=0.2, col="red")
            polygon(x,y,border="red", lwd=0.1)
          }
        }
      })
      ##Distance by species
      output$plot6 <- renderPlot({
        sp=paste(RES()$Complete$outsp)
        TAB<-trees2matrices.Distatis(trees, distance=choice())
        nam<-trees[[1]]$tip.label
        #TAB<-lapply(TAB, function(x,y) x[y,y],y=nam)
        par(mfrow=c(ceiling(length(trees[[1]]$tip.label)/5),5))
        par(mar=c(0,0,0,0))
        par(oma=c(0,0,0,0))
        for (j in 1:length(nam)) { ##for each species
          SP<-nam[j]
          T1<-lapply(TAB, function(x) (x[SP,nam]))
          T1m<-matrix(unlist(T1), nrow=length(trees), byrow=TRUE)
          ##T1m gives 1 plot corresponding to "Kla" for each gene.
          Means.T1m<-apply(T1m, 2, mean)
          ##CIRCLE:
          xc<-rep(1,length(nam)+1)*cos(seq(0,2*pi,length.out=length(nam)+1))
          yc<-rep(1,length(nam)+1)*sin(seq(0,2*pi,length.out=length(nam)+1))
          ##we check angles
          alphas<-seq(0,2*pi,length.out=length(nam)+1)
          alphas<-alphas[1:length(nam)]
          xc<-xc[1:length(nam)]
          yc<-yc[1:length(nam)]
          alphas<-alphas[1:length(nam)]
          ##for each gene, the ray is given by the proportion:
          GENEi<-NULL
          plot(4*xc,4*yc,type="n", xlim=c(-4,4), ylim=c(-4,4), frame.plot=FALSE, axes=FALSE, xlab="", ylab="")
          if (SP%in%sp){
            colo = "red"
          }
          else{
            colo= "black"
          }
          text(4*xc,4*yc,labels=nam, col="light grey", col = colo)
          for (i in 1:length(trees)) {
            genei<-T1m[i,]/Means.T1m
            genei[is.na(genei)]<-1
            GENEi<-c(GENEi, genei)
            x<-genei*cos(alphas)
            y<-genei*sin(alphas)
            x[is.na(x)]<-0
            y[is.na(y)]<-0
            polygon(xc,yc, border="light grey", lwd=0.54)
            polygon(x,y,border="red", lwd=0.8)
            text(-3.5,-3.5,SP,cex=2)
          }
        }
      })
    }
      
      #################################################################################################event = DETECTION OUTLIER OFF
    else{
      output$outputPhylter1 <- renderText ({
        ""
      })
        output$outputPhylter2 <- renderText ({
        ""
      })
      output$outputPhylter3 <- renderText ({
        ""
      })
     
      output$plot3 <- renderPlot({
        participant.colors = "grey50"
        plot(G, col = participant.colors, cex = 3.5)
        text(G, labels = rownames(G))
      })
      ##Distance by species
      output$plot6 <- renderPlot({
        TAB<-trees2matrices.Distatis(trees, distance=choice())
        nam<-trees[[1]]$tip.label
        #TAB<-lapply(TAB, function(x,y) x[y,y],y=nam)
        par(mfrow=c(ceiling(length(trees[[1]]$tip.label)/5),5))
        par(mar=c(0,0,0,0))
        par(oma=c(0,0,0,0))
        for (j in 1:length(nam)) { ##for each species
          SP<-nam[j]
          T1<-lapply(TAB, function(x) (x[SP,nam]))
          T1m<-matrix(unlist(T1), nrow=length(trees), byrow=TRUE)
          ##T1m gives 1 plot corresponding to "Kla" for each gene.
          Means.T1m<-apply(T1m, 2, mean)
          ##CIRCLE:
          xc<-rep(1,length(nam)+1)*cos(seq(0,2*pi,length.out=length(nam)+1))
          yc<-rep(1,length(nam)+1)*sin(seq(0,2*pi,length.out=length(nam)+1))
          ##we check angles
          alphas<-seq(0,2*pi,length.out=length(nam)+1)
          alphas<-alphas[1:length(nam)]
          xc<-xc[1:length(nam)]
          yc<-yc[1:length(nam)]
          alphas<-alphas[1:length(nam)]
          ##for each gene, the ray is given by the proportion:
          GENEi<-NULL
          plot(4*xc,4*yc,type="n", xlim=c(-4,4), ylim=c(-4,4), frame.plot=FALSE, axes=FALSE, xlab="", ylab="")
          text(4*xc,4*yc,labels=nam, col="light grey")
          for (i in 1:length(trees)) {
            genei<-T1m[i,]/Means.T1m
            genei[is.na(genei)]<-1
            GENEi<-c(GENEi, genei)
            x<-genei*cos(alphas)
            y<-genei*sin(alphas)
            x[is.na(x)]<-0
            y[is.na(y)]<-0
            polygon(xc,yc, border="light grey", lwd=0.54)
            polygon(x,y,border="red", lwd=0.8)
            text(-3.5,-3.5,SP,cex=2)
          }
        }
      })
      ##Distance by genes
      output$plot7 <- renderPlot({
        Geneslist = input$Geneslist
        Geneslist = strsplit(Geneslist, "-")
        TAB<-trees2matrices.Distatis(trees, distance=choice())
        nam<-trees[[1]]$tip.label
        #TAB<-lapply(TAB, function(x,y) x[y,y],y=nam)
        par(mfrow=c(5,4))
        par(mar=c(1,1,1,1))
        par(oma=c(0,0,0,0))
        i1 = as.integer(Geneslist[[1]][1])
        i2 = as.integer(Geneslist[[1]][2])
        for (i in i1:i2) {
          plot(0,0,type="n", xlim=c(-4,4), ylim=c(-4,4), frame.plot=FALSE, axes=FALSE, xlab="", ylab="", main=i, col.main="black", cex.main = 1.5)
          for (j in 1:length(nam)) { ##for each speciew
            SP<-nam[j]
            T1<-lapply(TAB, function(x,y) (x[SP,nam]))
            T1m<-matrix(unlist(T1), nrow=length(trees), byrow=TRUE)
            Means.T1m<-apply(T1m, 2, mean)
            genei<-T1m[i,]/Means.T1m
            xc<-rep(1,length(nam)+1)*cos(seq(0,2*pi,length.out=length(nam)+1))
            yc<-rep(1,length(nam)+1)*sin(seq(0,2*pi,length.out=length(nam)+1))
            ##we check angles
            alphas<-seq(0,2*pi,length.out=length(nam)+1)
            alphas<-alphas[1:length(nam)]
            x<-genei*cos(alphas)
            y<-genei*sin(alphas)
            x[is.na(x)]<-0
            y[is.na(y)]<-0
            polygon(xc,yc, border="light grey", lwd=0.5)
            points(x,y,pch=19, cex=0.2, col="red")
            polygon(x,y,border="red", lwd=0.1)
          }
        }
      })  
    }
  })
  ##Species

  observe({
    updateSelectInput(session, "selectSpecies", selected= "all", choices = append(append(outVar(),"mean",after=0),"all", after=0))
  })
  observeEvent(input$selectSpecies, {
    S=input$selectSpecies
    if (S == "all"){
      output$plot2 <- renderPlot({
        par(mfrow=c(ceiling(length(trees[[1]]$tip.label)/4),4))
        par(mar=c(2,2,2,2))
        part.design <- diag(dim(PF)[3])
        participant.colors <- as.matrix(createColorVectorsByDesign(part.design)$oc)
        for (sp in 1:length(trees[[1]]$tip.label)){
          to.plot <- t(PF[sp, , ])
          center.point <- FS[sp, c(1, 2)]
          center.rep <- matrix(center.point, dim(PF)[3], 2, byrow = TRUE)
          bound.mat <- rbind(center.rep, to.plot[, c(1, 2)])
          bound.mat <- bound.mat[as.vector(t(matrix(seq(1, nrow(bound.mat)), ncol = 2))), ]
          plot(to.plot, main = dimnames(PF)[[1]][sp], cex.main= 1.5, col = participant.colors, cex=3.2)
          points(bound.mat, type = "l", lty = 1, lwd = 1, col = "grey70")
          text(to.plot, labels = rownames(t(PF[sp, , ])))
        }
      })
    }
    if(S == "mean"){
      output$plot2 <- renderPlot({
        item.design <- diag(dim(FS)[1])
        item.colors <- as.matrix(createColorVectorsByDesign(item.design)$oc)
        real.minimum <- min(FS)
        real.maximum <- max(FS)
        real.value <- max(c(abs(real.minimum), abs(real.maximum)))
        plot(FS, col = item.colors,axes=F, cex = 3.5, ylim=c(-real.value,real.value),xlim=c(-real.value,real.value))
        axis(1, pos=0,lwd.ticks=0.5)
        axis(2, pos=0,lwd.ticks=0.5)
        text(FS, rownames(FS))
      })
    }
    if(S != "mean" && S != "all"){ 
      output$plot2 <- renderPlot({
        part.design <- diag(dim(PF)[3])
        participant.colors <- as.matrix(createColorVectorsByDesign(part.design)$oc)
        to.plot <- t(PF[S, , ])
        center.point <- FS[S, c(1, 2)]
        center.rep <- matrix(center.point, dim(PF)[3], 2, byrow = TRUE)
        bound.mat <- rbind(center.rep, to.plot[, c(1, 2)])
        bound.mat <- bound.mat[as.vector(t(matrix(seq(1, nrow(bound.mat)), ncol = 2))), ]
        plot(to.plot, main = dimnames(PF)[[1]][S], cex.main= 1.5, col = participant.colors, cex=3.2)
        points(bound.mat, type = "l", lty = 1, lwd = 1, col = "grey70")
        text(to.plot, labels = rownames(t(PF[S, , ])))
      })
    }
  })
  observeEvent(input$BackClick2, {
    updateSelectInput(session, "selectSpecies", selected= "all")
  })

  ##Gene clusters
    output$plot5 <- renderPlot({
      plot(hclust(as.dist(C)))
    })
  ##WR
    output$plot4 <- renderPlot({
      WR<-WR()
      pl = plot2WR(WR)
      pl + coord_cartesian(xlim = ranges$x, ylim = ranges$y, expand = FALSE)
      #pl + coord_fixed(ratio=1/5)
    })
    observeEvent(input$plot4dbclick, {
      brush <- input$plot4_brush
      if (!is.null(brush)) {
        ranges$x <- c(brush$xmin, brush$xmax)
        ranges$y <- c(brush$ymin, brush$ymax)
        
      } else {
        ranges$x <- NULL
        ranges$y <- NULL
      }
    })
  })
}

