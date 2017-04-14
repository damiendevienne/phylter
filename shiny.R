library(shiny)
library(shinythemes)
#library(plotrix)

source("/home/aurore/Documents/Phylter/pmcoa.R")
source("/home/aurore/Documents/Phylter/PhylteR.R")

#source("F:/stageLBBE/Phylter-R/pmcoa.R")
#source("F:/stageLBBE/Phylter-R/PhylteR.R")

#trees = read.tree("/home/aurore/Documents/Phylter/viz/Aguileta-et-al-2008_TREES.txt", keep.multi = TRUE)
#trees=read.tree(file="/home/aurore/Documents/Phylter/trees/rose/test.phy", keep.multi = TRUE)
ui<-fluidPage(
  theme = shinytheme("flatly"),
  sidebarLayout(
    sidebarPanel(
      wellPanel(
        titlePanel("Data input"),
        fileInput(inputId = "trees", label = "add a tree file"),
        radioButtons(inputId="choice",label = "Method", choiceNames = c("nodal","patristic"), choiceValues = c("nodal","patristic"))
      ),
      wellPanel(
        titlePanel("PhylteR"),
        radioButtons(inputId="Outlier",label = "Detection", choiceNames = c("off","on"), choiceValues = c("off","on"), selected = "off"),
        sliderInput(inputId="k", label = "select k", value = 3, min=1, max=6,round=FALSE,step=0.1)
        ),
      wellPanel(
        titlePanel("Visualisation"),
        actionButton(inputId="WR",label = "WR"),
        actionButton(inputId="species",label = "SpeciesByGenes"),
        actionButton(inputId="species2",label = "Species"),
        actionButton(inputId="genesG",label = "Genes"),
        actionButton(inputId="genes2",label = "Genes2"),
        actionButton(inputId="genesCortrees",label = "Genes Clusters"),
        selectInput(inputId = "selectSpecies",label = "see a particular specie","")
      )
    ),
    mainPanel(
      textOutput(outputId= "outputPhylter1"),
      textOutput(outputId= "outputPhylter2"),
      textOutput(outputId= "outputPhylter3"),
      plotOutput(outputId = "plot1", height = "1000px",dblclick = "BackClick")
    )
  )
)

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
  
  WR <- reactive({
    Dist2WR(dist())
  })
  trees <- reactive({
    read.tree(input$trees$datapath, keep.multi = TRUE)
  })
  
###########-----------------------------------

  ##event = DETECTION OUTLIER ON
  
  observeEvent(outliers(), {
    if (outliers()== "on"){
    trees = trees()
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
    observeEvent(input$genesG, {
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
      output$plot1 <- renderPlot({
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
    })
    }
    else{
      ##event = DETECTION OUTLIER OFF
        output$outputPhylter1 <- renderText ({
          ""
        })
        output$outputPhylter2 <- renderText ({
          ""
        })
        output$outputPhylter3 <- renderText ({
          ""
        })
        observeEvent(input$genesG, {
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
          output$plot1 <- renderPlot({
            participant.colors = "grey50"
            plot(G, col = participant.colors, cex = 3.5)
            text(G, labels = rownames(G))
          })
        })
    }
  })
  
  observeEvent(input$genesCortrees, {
    C = dist()$res4Cmat$C
    output$plot1 <- renderPlot({
      plot(hclust(as.dist(C)))
    })
  })
  observeEvent(input$species, {
    trees = trees()
    FS=dist()$res4Splus$F
    PF=dist()$res4Splus$PartialF
    if(is.null(dimnames(PF)[[3]])){
      if(!is.null(names(trees))){
        dimnames(PF)[[3]] = names(trees)
      }
      else {
        dimnames(PF)[[3]] = as.character(c(1:length(trees)))
      }
    }
    observe({
      updateSelectInput(session, "selectSpecies", selected= "all", choices = append(append(outVar(),"mean",after=0),"all",after=0))
    })
    observeEvent(input$selectSpecies, {
      trees = trees()
        S=input$selectSpecies
       if (S == "all"){
          output$plot1 <- renderPlot({
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
        output$plot1 <- renderPlot({
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
          output$plot1 <- renderPlot({
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
      observeEvent(input$BackClick, {
        updateSelectInput(session, "selectSpecies", selected= "all")
      })
    })

  observeEvent(input$WR, {
    output$plot1 <- renderPlot({
      #plot.2WR(t(WR))
      par(mar=c(1,1,1,1))
      WR<-normalize(WR())
      testlat<-require(lattice)
      ok<-levelplot(WR, xlab="Species", ylab="Genes", col.regions=grey(seq(1,0,length.out=100)), scales = list(x = list(rot = 90)))
      plot(ok)
      
    })
  })
  observeEvent(input$Species2, { 
    output$plot1 <- renderPlot({
      trees = trees()
      nam = trees[[1]]$tip.label
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
        for (i in 1:length(tr)) {
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
    })
  })
}


shinyApp(ui=ui,server =server)
