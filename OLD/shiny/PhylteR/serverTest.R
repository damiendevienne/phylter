#server options : 5 Mo maximum.
options(shiny.maxRequestSize = 4*1024^2)

server <-function(input,output,session){

 ##reactives
  outVar <- reactive({
    TAB=mat()
    sps<-rownames(TAB[[1]])
    if (!is.null(sps)){
      vars <- all.vars(parse(text = sps))
    }
  })
  outliers <- reactive({
    input$Outlier
  })
  RES <- reactive({
    withProgress(message = 'PhylteR is running', value = 0, {
    trees = trees()
      if(is.null(names(trees))){
        names(trees)=as.character(c(1:length(trees)))
      }
      incProgress(1)
      PhylteR(trees, distance = choice(), k=k(), thres = 0.5)
    })
  })
  outgn <- reactive({
    paste(RES()$Complete$outgn)
  })
  outsp <- reactive({
    paste(RES()$Complete$outsp)
  })
  outcell <- reactive({
    if(!is.null(RES()$CellByCell$outcell)){
      if(!is.null(nrow(RES()$CellByCell$outcell))){
        data.frame(paste(RES()$CellByCell$outcell[,1]),paste(RES()$CellByCell$outcell[,2]))
      }
      #gestion de l'exeption qui fait que quand on a qu'un seul outlierCell la ligne unique n'est pas détecté comme une dataframe
      else {
        data.frame(paste(RES()$CellByCell$outcell[[1]]),paste(RES()$CellByCell$outcell[[2]]))
      }
    }
    else {
      NULL
    }
  })
  mat <- reactive({
    withProgress(message = 'imputation', value = 0, {
      incProgress(1)
      trees = trees()
      if(is.null(names(trees))){
        names(trees)=as.character(c(1:length(trees)))
      }
      mat <- trees2matrices(trees, choice())
      impPCA.multi(mat,maxiter=10)
    })
  })
  dist <- reactive({
    mat <- mat()
    mat2Dist(mat)
  })
  choice <- reactive({
    input$choice
  })
  k <- reactive({
    input$k
  })
  TreesList <- reactive({
    input$Treeslist
  })
  WR <- reactive({
    Dist2WR(dist())
  })
  trees <- reactive({
    withProgress(message = 'Trees are loading', value = 0, {
      incProgress(1)
      read.tree(input$trees$datapath, keep.multi = TRUE)
    })
  })
  PlotHeight = reactive(
    50*(length(outVar()))
  )
  ranges <- reactiveValues(x = NULL, y = NULL)
  #end of reactives

  #download file test
  output$downloadData <- downloadHandler(
    filename = "Aguileta-et-al-2008_TREES.txt",
    content = function(file) {
      file.copy("www/Aguileta-et-al-2008_TREES.txt", file)
    },
    contentType = "txt"
  )

  #What happens when we add datas
  observeEvent(input$trees, {
    
    output$DataSetInfo <- renderText ({
      paste("This data set has ", length(trees()), " genes and ", length(length(outVar())), " species", sep="")
    })

    #activation of tabs
    session$sendCustomMessage('activeNavs', 'visualize some trees')
    session$sendCustomMessage('activeNavs', 'visualize species on distatis compromise')
    session$sendCustomMessage('activeNavs', 'visualize genes')
    session$sendCustomMessage('activeNavs', 'visualize distances')
    session$sendCustomMessage('activeNavs', 'visualize 2WR')
    
    #functions from PhylteR

    FS=dist()$res4Splus$F
    PF=dist()$res4Splus$PartialF
    C = dist()$res4Cmat$C
    if(is.null(dimnames(PF)[[3]])){
      if(!is.null(names(trees()))){
        dimnames(PF)[[3]] = names(trees())
      }
      else {
        dimnames(PF)[[3]] = as.character(c(1:length(trees())))
      }
    }
    G=dist()$res4Cmat$G
    if(is.null(rownames(G))){
      if(!is.null(names(trees()))){
        rownames(G)<-names(trees())
      }
      else{
        rownames(G)<-as.character(c(1:length(trees())))
      }
    }
    
    ##Update Genetrees. Genelist from the panel "visualize distances between genes". Genes are spilts in partition of 50 genes
    
    trees = trees()
    if(is.null(names(trees))){
      names(trees)=as.character(c(1:length(trees)))
    }
    n=names(trees)
    c=1
    l1 = 50
    Partition = list()
    while (l1 < length(n)){
      p=n[l1]
      k=paste(c, "-", p, sep="")
      Partition = append(Partition, k)
      l1 =l1+50
      c=c+50
    }
    if (l1 >= length(n)){
      k=paste(c, "-", length(n), sep="")
      Partition = append(Partition, k)
    }
    observe({
      updateSelectInput(session, "Geneslist", choices = Partition)
    })

#Several tabs and plots are differents if the detection is on or off
    observeEvent(outliers(), {
      
      ##event = DETECTION OUTLIER ON
      if (outliers()== "on"){
        trees = trees()
        if(is.null(names(trees))){
          names(trees)=as.character(c(1:length(trees)))
        }
        
        ##Visualize some trees tab. We can type several trees in the box "Sumit trees", all separated by a "-" and trees are ploted.
        observeEvent(input$SumitTrees,{
          List = TreesList()
          S = strsplit(List, ",")
          S = S[[1]]
          output$plot1 <- renderPlot(height = 500*length(S), {
            withProgress(message = 'Making plot', value = 0, {
              if (length(S)==1){par(mfrow = c(1,1))}
              if (length(S)>=2 && length(S)<=8){par(mfrow = c(ceiling(length(S)/2),2))}
              if (length(S)>8 && length(S)<=30){par(mfrow = c(ceiling(length(S)/5),5))}
              if (length(S)>30){par(mfrow = c(ceiling(length(S)/8),8))}
              par(mar = c(0,0,1,0))
              for (i in 1:length(S)){
                incProgress(1/length(S))
                gene.colors <- "black"
                species.colors <- length(outVar())
                if (length(outgn())!=0){
                  for (j in 1:length(outgn())){
                    if (outgn()[j]== S[i]){
                      gene.colors = "red"
                    }
                  }
                }
                if (length(outsp())!=0){
                  for (j in 1:length(outsp())){
                      species.colors[which(species.colors==outsp()[j])] = "darkgreen"
                  }
                }
                if (nrow(outcell())!=0){
                  for (j in 1:length(outcell()[,2])){
                    if(outcell()[j,2]==S[i]){
                      species.colors[which(species.colors==outcell()[j,1])] = "blue"
                    }
                  }
                }
                for (k in 1:length(species.colors)){
                  if(species.colors[k] != "darkgreen" && species.colors[k] != "blue"){
                    species.colors[k] <- "black"
                  }
                }
                ape::plot.phylo(trees[[S[i]]],  edge.color = gene.colors, tip.color = species.colors, main = S[i])
              }
            })
          })
        })
        ##Observe species:
        observe({
          updateSelectInput(session, "selectSpecies", selected= "all", choices = append(append(outVar(),"mean",after=0),"all", after=0))
        })
        observeEvent(input$selectSpecies, {
          S=input$selectSpecies
          #all species in different plot. All in the same scale.
          if (S == "all"){
            output$plot2 <- renderPlot(height = PlotHeight(), {
              withProgress(message = 'Making plot', value = 0, {
                par(mfrow=c(ceiling(length(outVar())/4),4))
                par(mar=c(2,2,2,2))
                part.design <- diag(dim(PF)[3])

                for (sp in 1:length(outVar())){
                  
                  participant.colors <- rownames(G)
                  fill.colors <- rownames(G)
                  if (length(outgn())!=0){
                    for (j in 1:length(outgn())){
                      participant.colors[which(participant.colors==outgn()[j])] = "red1"
                      fill.colors[which(fill.colors==outgn()[j])] = "pink"
                    }
                  }
                  if (length(outcell()[,1])!=0){
                    for (j in 1:length(outcell()[,2])){
                      if(outcell()[j,1]==outVar()[sp]){
                        participant.colors[which(participant.colors==outcell()[j,2])] = "blue"
                        fill.colors[which(fill.colors==outcell()[j,2])] = "paleturquoise"
                      }
                    }
                  }
                  for (color in 1:length(participant.colors)){
                    if(participant.colors[color] != "red1" && participant.colors[color] != "blue"){
                      participant.colors[color] <- "grey50" 
                      fill.colors[color]="white"
                    }
                  }
                  
                  incProgress(1/length(outVar()))
                  to.plot <- t(PF[sp, , ])
                  center.point <- FS[sp, c(1, 2)]
                  center.rep <- matrix(center.point, dim(PF)[3], 2, byrow = TRUE)
                  bound.mat <- rbind(center.rep, to.plot[, c(1, 2)])
                  bound.mat <- bound.mat[as.vector(t(matrix(seq(1, nrow(bound.mat)), ncol = 2))), ]
                  coltitle="black"
                  if(length(outsp())!=0){
                    for (out in 1:length(outsp())){
                      if (outsp()[out] == outVar()[sp]){
                        coltitle = "darkgreen"
                      }
                    }
                  }
                  plot(to.plot, main = dimnames(PF)[[1]][sp], col.main=coltitle, cex.main= 1.5, col = participant.colors, cex=3.2, xlim = range(t(PF[, 1, ])), ylim = range(t(PF[, 2, ])), pch = 21, bg=fill.colors)
                  points(bound.mat, type = "l", lty = 1, lwd = 1, col = "grey70")
                  text(to.plot, labels = rownames(t(PF[sp, , ])))
                }
              })
            })
          }
          #mean position of each specie on a graph
          if(S == "mean"){
            output$plot2 <- renderPlot({
              withProgress(message = 'Making plot', value = 0, {
                incProgress(0)
                item.design <- diag(dim(FS)[1])

                participant.colors <- rownames(FS)
                fill.colors <- rownames(FS)
                if (length(outsp())!=0){
                  for (j in 1:length(outsp())){
                    participant.colors[which(participant.colors==outsp())[j]] = "darkgreen"
                    fill.colors[which(fill.colors==outsp()[j])] = "palegreen"
                  }
                  participant.colors[which(participant.colors!="darkgreen")] = "grey50"
                  fill.colors[which(fill.colors!="palegreen")] = "white"
                }
                else{
                  participant.colors = "grey50"
                  fill.colors="white"
                }
                
                real.minimum <- min(FS)
                real.maximum <- max(FS)
                real.value <- max(c(abs(real.minimum), abs(real.maximum)))
                plot(FS, col = participant.colors, axes=F, cex = 3.5, ylim=c(-real.value,real.value),xlim=c(-real.value,real.value),pch = 21,bg=fill.colors)
                axis(1, pos=0,lwd.ticks=0.5)
                axis(2, pos=0,lwd.ticks=0.5)
                text(FS, rownames(FS))
              })
            })
          }
          #single specie selected
          if(S != "mean" && S != "all"){
            output$plot2 <- renderPlot({
              withProgress(message = 'Making plot', value = 0, {
                incProgress(0)
                part.design <- diag(dim(PF)[3])

                participant.colors <- rownames(G)
                fill.colors <- rownames(G)
                
                if (length(outgn())!=0){
                  for (j in 1:length(outgn())){
                    participant.colors[which(participant.colors==outgn()[j])] = "red1"
                    fill.colors[which(fill.colors==outgn()[j])] = "pink"
                  }
                }
                  if (length(outcell()[,2])!=0){
                    for (j in 1:length(outcell()[,2])){
                        if(outcell()[j,1]==S){
                          participant.colors[which(participant.colors==outcell()[j,2])] = "blue"
                          fill.colors[which(fill.colors==outcell()[j,2])] = "paleturquoise"
                      }
                    }
                  }
                  for (color in 1:length(participant.colors)){
                    if(participant.colors[color] != "red1" && participant.colors[color] != "blue"){
                      participant.colors[color] <- "grey50" 
                      fill.colors[color]="white"
                    }
                  }
                to.plot <- t(PF[S, , ])
                center.point <- FS[S, c(1, 2)]
                center.rep <- matrix(center.point, dim(PF)[3], 2, byrow = TRUE)
                bound.mat <- rbind(center.rep, to.plot[, c(1, 2)])
                bound.mat <- bound.mat[as.vector(t(matrix(seq(1, nrow(bound.mat)), ncol = 2))), ]
                plot(to.plot, main = dimnames(PF)[[1]][S], cex.main= 1.5, col = participant.colors, cex=3.2,pch = 21,bg=fill.colors)
                points(bound.mat, type = "l", lty = 1, lwd = 1, col = "grey70")
                text(to.plot, labels = rownames(t(PF[S, , ])))
              })
            })
          }
        })
        #Back clicking = return on the "all" plot.
        observeEvent(input$BackClick2, {
          updateSelectInput(session, "selectSpecies", selected= "all")
        })
        
        #Update the text zone with outliers detected:
        #outgn
        output$outputPhylter1 <- renderText ({
          if(length(outgn())!=0){
            gn=paste(outgn(), ";")
            append (gn, "outlier genes: ", after=0)
          }
          else{
            "no outlier gene detected"
          }
        })
        #sp
        output$outputPhylter2 <- renderText ({
          if(length(outsp())!=0){
            sp=paste(outsp(), ";")
            append (sp, "outlier species: ", after=0)
          }
          else{
            "no outlier species detected"
          }

        })
        #cell
        output$outputPhylter3 <- renderText ({
          if(nrow(outcell())!=0){
            ce = paste(outcell()[,1],outcell()[,2], ";")
            append (ce, "outlier gene/species: ", after=0)
          }
          else{
            "no outlier cell detected"
          }
        })
    #visualize genes: genes are red if outliers, grey if not
      output$plot3 <- renderPlot({
        withProgress(message = 'Making plot', value = 0, {
        participant.colors <- rownames(G)
        incProgress(1)
        if (length(outgn())!=0){
          for (j in 1:length(outgn())){
            participant.colors[which(participant.colors==outgn()[j])] = "red1"
          }
          participant.colors[which(participant.colors!="red1")] = "grey50"
        }
        else{
          participant.colors = "grey50"
        }
        plot(G, col = participant.colors, cex = 3.5)
        text(G, labels = rownames(G))
        })
      })
      ##visualize distances between genes : outliers detected are red.
      output$plot7 <- renderPlot({
        withProgress(message = 'Making plot', value = 0, {
        Geneslist = input$Geneslist
        Geneslist = strsplit(Geneslist, "-")
        TAB<-mat()
        nam<-rownames(TAB[[1]])
        #TAB<-lapply(TAB, function(x,y) x[y,y],y=nam)
        par(mfrow=c(5,10))
        par(mar=c(0,0,0,0))
        par(oma=c(0,0,0,0))
        i1 = as.integer(Geneslist[[1]][1])
        i2 = as.integer(Geneslist[[1]][2])
        
        listx=vector()
        listy=vector()
        for (i in i1:i2) {
          for (j in 1:length(nam)) { ##for each speciew
            SP<-nam[j]
            T1<-lapply(TAB, function(x,y) (x[SP,nam]))
            T1m<-matrix(unlist(T1), nrow=length(trees), byrow=TRUE)
            Means.T1m<-apply(T1m, 2, mean)
            genei<-T1m[i,]/Means.T1m
            alphas<-seq(0,2*pi,length.out=length(nam)+1)
            alphas<-alphas[1:length(nam)]
            x<-genei*cos(alphas)
            y<-genei*sin(alphas)
            x[is.na(x)]<-0
            y[is.na(y)]<-0
            listx = append(listx,x)
            listy = append(listy,y)
          }
        }
        
        for (i in i1:i2) {
          incProgress(1/length(Geneslist))
          if (i%in%outgn()){
            colo = "red"
          }
          else{
            colo = "black"
          }
          plot(0,0,type="n", xlim=c(min(listx),max(listx)), ylim=c(min(listy),max(listy)), frame.plot=FALSE, axes=FALSE, xlab="", ylab="", col.main=colo, cex.main = 1.5)
          title(i,line=-5)
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
    })
      ##Visualize distance between species
      output$plot6 <- renderPlot(height = 50*length(outVar()),{
        withProgress(message = 'Making plot', value = 0, {
        sp=paste(RES()$Complete$outsp)
        TAB<-mat()
        nam<-rownames(TAB[[1]])
        #TAB<-lapply(TAB, function(x,y) x[y,y],y=nam)
        par(mfrow=c(ceiling(length(outVar())/5),5))
        par(mar=c(0,0,0,0))
        par(oma=c(0,0,0,0))
        
        listx =vector()
        listy =vector()
        for (j in 1:length(nam)) {
          GENEi<-NULL
          SP<-nam[j]
          T1<-lapply(TAB, function(x) (x[SP,nam]))
          T1m<-matrix(unlist(T1), nrow=length(trees), byrow=TRUE)
          ##T1m gives 1 plot corresponding to "Kla" for each gene.
          Means.T1m<-apply(T1m, 2, mean)
          ##we check angles
          alphas<-seq(0,2*pi,length.out=length(nam)+1)
          alphas<-alphas[1:length(nam)]
          alphas<-alphas[1:length(nam)]
          for (i in 1:length(trees)) {
            genei<-T1m[i,]/Means.T1m
            genei[is.na(genei)]<-1
            GENEi<-c(GENEi, genei)
            x<-genei*cos(alphas)
            y<-genei*sin(alphas)
            x[is.na(x)]<-0
            y[is.na(y)]<-0
            listx = append(listx,x)
            listy = append(listy,y)
          }
        }
        
        for (j in 1:length(nam)) { ##for each species
          incProgress(1/length(nam))
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
          listx=vector()
          listy=vector()
          for (i in 1:length(trees)) {
            genei<-T1m[i,]/Means.T1m
            genei[is.na(genei)]<-1
            GENEi<-c(GENEi, genei)
            x<-genei*cos(alphas)
            y<-genei*sin(alphas)
            x[is.na(x)]<-0
            y[is.na(y)]<-0
            listx=append(listx,x)
            listy=append(listy,y)
          }
          plot(4*xc,4*yc,type="n", xlim=c(min(listx),max(listx)), ylim=c(min(listy),max(listy)), frame.plot=FALSE, axes=FALSE, xlab="", ylab="")
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
            if (SP%in%sp){
              colo = "red"
            }
            else{
              colo= "black"
            }
            text(-3.5,-3.5,SP,cex=2, col = colo)
          }
        }
      })
    })
      ##Visualize the 2WR matrix
      output$plot4 <- renderPlot({
        withProgress(message = 'Making plot', value = 0, {
          incProgress(1)
          
          mat = mat()
          
          genes = names(trees)
          species = rownames(mat[[1]])
          genesDF = vector()
          for (i in 1:length(species)){
            genesDF = append(genesDF,genes)
          }
          speciesDF = vector()
          for (i in 1:length(genes)){
            speciesDF = append(speciesDF,species)
          }
          dataoutlier = data.frame(
            gn = sort(as.numeric(genesDF),decreasing = FALSE),
            sp = speciesDF
          )
          
          ##outgn
          gn=NULL
          sp=NULL
          gn = cbind(rep(as.vector(outgn()),each=length(species)))
          sp = rep(species,length(outgn()))
          datag = data.frame(gn,sp)
          
          ##outsp
          gn=NULL
          sp=NULL
          gn = cbind(rep(as.vector(outsp()),each=length(genes)))
          sp = rep(genes,length(outsp()))
          datas = data.frame(gn,sp)
          
          ##outcell
          gn=NULL
          sp=NULL
          gn = as.vector(outcell()[,2])
          sp = as.vector(outcell()[,1])
          
          datac = data.frame(gn,sp)
          
          dataoutlier = rbind(datag, datas, datac)
          dataoutlier=data.frame(dataoutlier,y=1)
          
          WR<-WR()
          pl = plot2WR(WR)
          #pl = pl + ggplot(aes(data=dataoutlier,x=gn,y=sp,color=factor(y)))
          pl + coord_cartesian(xlim = ranges$x, ylim = ranges$y, expand = FALSE)
          
        })
      })
    }

      ##event = DETECTION OUTLIER OFF
    else{
      trees = trees()
      if(is.null(names(trees))){
        names(trees)=as.character(c(1:length(trees)))
      }
      
      ##Visualize some trees tab. We can type several trees in the box "Sumit trees", all separated by a "-" and trees are ploted.
      observeEvent(input$SumitTrees,{
        List = TreesList()
        S = strsplit(List, ",")
        S = S[[1]]
        output$plot1 <- renderPlot(height = 500*length(S),{
          withProgress(message = 'Making plot', value = 0, {
            if (length(S)==1){par(mfrow = c(1,1))}
            if (length(S)>=2 && length(S)<=8){par(mfrow = c(ceiling(length(S)/2),2))}
            if (length(S)>8 && length(S)<=30){par(mfrow = c(ceiling(length(S)/5),5))}
            if (length(S)>30){par(mfrow = c(ceiling(length(S)/8),8))}
            par(mar = c(0,0,1,0))
            for (i in 1:length(S)){
              incProgress(1/length(S))
              ape::plot.phylo(trees[[S[i]]], main = S[i])
            }
          })
        })
      })
      
      ##Observe species:
      observe({
        updateSelectInput(session, "selectSpecies", selected= "all", choices = append(append(outVar(),"mean",after=0),"all", after=0))
      })
      observeEvent(input$selectSpecies, {
        S=input$selectSpecies
        #all species in different plot. All in the same scale.
        if (S == "all"){
          output$plot2 <- renderPlot(height = PlotHeight(),{
            withProgress(message = 'Making plot', value = 0, {
              par(mfrow=c(ceiling(length(outVar())/4),4))
              par(mar=c(2,2,2,2))
              part.design <- diag(dim(PF)[3])
              participant.colors <- "grey70"
              for (sp in 1:length(outVar())){
                incProgress(1/length(outVar()))
                to.plot <- t(PF[sp, , ])
                center.point <- FS[sp, c(1, 2)]
                center.rep <- matrix(center.point, dim(PF)[3], 2, byrow = TRUE)
                bound.mat <- rbind(center.rep, to.plot[, c(1, 2)])
                bound.mat <- bound.mat[as.vector(t(matrix(seq(1, nrow(bound.mat)), ncol = 2))), ]
                plot(to.plot, main = dimnames(PF)[[1]][sp], cex.main= 1.5, col = participant.colors, cex=3.2, xlim = range(t(PF[, 1, ])), ylim = range(t(PF[, 2, ])))
                points(bound.mat, type = "l", lty = 1, lwd = 1, col = "grey70")
                text(to.plot, labels = rownames(t(PF[sp, , ])))
              }
            })
          })
        }
        #mean position of each specie on a graph
        if(S == "mean"){
          output$plot2 <- renderPlot({
            withProgress(message = 'Making plot', value = 0, {
              incProgress(0)
              item.design <- diag(dim(FS)[1])
              item.colors <- "grey70"
              real.minimum <- min(FS)
              real.maximum <- max(FS)
              real.value <- max(c(abs(real.minimum), abs(real.maximum)))
              plot(FS, col = item.colors,axes=F, cex = 3.5, ylim=c(-real.value,real.value),xlim=c(-real.value,real.value))
              axis(1, pos=0,lwd.ticks=0.5)
              axis(2, pos=0,lwd.ticks=0.5)
              text(FS, rownames(FS))
            })
          })
        }
        #single specie selected
        if(S != "mean" && S != "all"){
          output$plot2 <- renderPlot({
            withProgress(message = 'Making plot', value = 0, {
              incProgress(0)
              part.design <- diag(dim(PF)[3])
              #participant.colors <- as.matrix(createColorVectorsByDesign(part.design)$oc)
              participant.colors <- "grey70"
              to.plot <- t(PF[S, , ])
              center.point <- FS[S, c(1, 2)]
              center.rep <- matrix(center.point, dim(PF)[3], 2, byrow = TRUE)
              bound.mat <- rbind(center.rep, to.plot[, c(1, 2)])
              bound.mat <- bound.mat[as.vector(t(matrix(seq(1, nrow(bound.mat)), ncol = 2))), ]
              plot(to.plot, main = dimnames(PF)[[1]][S], cex.main= 1.5, col = participant.colors, cex=3.2)
              points(bound.mat, type = "l", lty = 1, lwd = 1, col = "grey70")
              text(to.plot, labels = rownames(t(PF[S, , ])))
            })
          })
        }
      })
      #Back clicking = return on the "all" plot.
      observeEvent(input$BackClick2, {
        updateSelectInput(session, "selectSpecies", selected= "all")
      })
  
      
      #update outliers text to blank
      output$outputPhylter1 <- renderText ({
        ""
      })
        output$outputPhylter2 <- renderText ({
        ""
      })
      output$outputPhylter3 <- renderText ({
        ""
      })
      #visualize genes : all genes are greys
      output$plot3 <- renderPlot({
        withProgress(message = 'Making plot', value = 0, {
        incProgress(1)
        participant.colors = "grey50"
        plot(G, col = participant.colors, cex = 3.5)
        text(G, labels = rownames(G))
        })
      })

      observe({        
        output$plot6 <- renderUI({
          TAB<-mat()
          nam<-rownames(TAB[[1]])
          plot_output_list <- lapply(1:length(nam), function(i) {
            plotname <- paste("plotS", i, sep="")
            plotOutput(plotname, height = 300, width = 500)
          })
          do.call(tagList, plot_output_list)
        })

     
      ###VISUALIZE SPECIES 
      
      TAB<-mat()
      nam<-rownames(TAB[[1]])
      #TAB<-lapply(TAB, function(x,y) x[y,y],y=nam)
      listx = vector()
      listy = vector()
      max_plots <- length(nam)
      for (j in 1:length(nam)) {
        GENEi<-NULL
        SP<-nam[j]
        T1<-lapply(TAB, function(x) (x[SP,nam]))
        T1m<-matrix(unlist(T1), nrow=length(trees), byrow=TRUE)
        Means.T1m<-apply(T1m, 2, mean)
        alphas<-seq(0,2*pi,length.out=length(nam)+1)
        alphas<-alphas[1:length(nam)]
        for (i in 1:length(trees)) {
          genei<-T1m[i,]/Means.T1m
          genei[is.na(genei)]<-1
          GENEi<-c(GENEi, genei)
          x<-genei*cos(alphas)
          y<-genei*sin(alphas)
          x[is.na(x)]<-0
          y[is.na(y)]<-0
          listx = append(listx,x)
          listy = append(listy,y)
        }
      }
      plot_list = vector()
        for (j in 1:max_plots) {
          local({
            my_i <- j
            plotname <- paste("plotS", my_i, sep="")
            SP<-nam[j]
            GENEi<-NULL
            T1<-lapply(TAB, function(x) (x[SP,nam]))
            T1m<-matrix(unlist(T1), nrow=length(trees), byrow=TRUE)
            ##T1m gives 1 plot corresponding to "Kla" for each gene.
            Means.T1m<-apply(T1m, 2, mean)
            ##we check angles
            alphas<-seq(0,2*pi,length.out=length(nam)+1)
            alphas<-alphas[1:length(nam)]
            ##CIRCLE:
            xc<-rep(1,length(nam)+1)*cos(seq(0,2*pi,length.out=length(nam)+1))
            yc<-rep(1,length(nam)+1)*sin(seq(0,2*pi,length.out=length(nam)+1))
            ##we check angles
            xc<-xc[1:length(nam)]
            yc<-yc[1:length(nam)]
            output[[plotname]] <- renderPlot({
                ##for each gene, the ray is given by the proportion:
                pt=plot(4*xc,4*yc,type="n", xlim=c(min(listx),max(listx)),ylim=c(min(listy),max(listy)), frame.plot=FALSE, axes=FALSE, xlab="", ylab="")
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
                plot_list=append(plot_list,pt)
              })  
            })
       grid.arrange(grobs=plot_list,ncol=10)
       }
      })
      
      ##Visialize distance by species
      output$plot6 <- renderPlot(height = 50*length(outVar()), {

      })
      ####Visialize distance by genes
      output$plot7 <- renderPlot({
        withProgress(message = 'Making plot', value = 0, {
        Geneslist = input$Geneslist
        Geneslist = strsplit(Geneslist, "-")
        TAB<-mat()
        nam<-rownames(TAB[[1]])
        #TAB<-lapply(TAB, function(x,y) x[y,y],y=nam)
        par(mfrow=c(5,10))
        par(mar=c(0,0,0,0))
        par(oma=c(0,0,0,0))
        i1 = as.integer(Geneslist[[1]][1])
        i2 = as.integer(Geneslist[[1]][2])
        
        listx=vector()
        listy=vector()
        for (i in i1:i2) {
          for (j in 1:length(nam)) { ##for each speciew
            SP<-nam[j]
            T1<-lapply(TAB, function(x,y) (x[SP,nam]))
            T1m<-matrix(unlist(T1), nrow=length(trees), byrow=TRUE)
            Means.T1m<-apply(T1m, 2, mean)
            genei<-T1m[i,]/Means.T1m
            alphas<-seq(0,2*pi,length.out=length(nam)+1)
            alphas<-alphas[1:length(nam)]
            x<-genei*cos(alphas)
            y<-genei*sin(alphas)
            x[is.na(x)]<-0
            y[is.na(y)]<-0
            listx = append(listx,x)
            listy = append(listy,y)
          }
        }
        
        for (i in i1:i2) {
          incProgress(1/length(Geneslist))
          plot(0,0,type="n", xlim=c(min(listx),max(listx)), ylim=c(min(listy),max(listy)), frame.plot=FALSE, axes=FALSE, xlab="", ylab="", col.main="black", cex.main = 1.5)
          title(i,line = -5)
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
      })
        ##Visualize the 2WR matrix
        output$plot4 <- renderPlot({
          withProgress(message = 'Making plot', value = 0, {
            incProgress(1)
            WR<-WR()
            pl = plot2WR(WR)
            pl + coord_cartesian(xlim = ranges$x, ylim = ranges$y, expand = FALSE)
        })
      })
    }
  })

#Evenements indépendant de la détection d'outliers  
    #zoom on 2WR matrix
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

