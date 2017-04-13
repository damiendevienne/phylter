library(shiny)
library(shinythemes)
library(shinyHeatmaply)
library(plotrix)
source("/home/aurore/Documents/Phylter/pmcoa.R")
source("/home/aurore/Documents/Phylter/PhylteR.R")

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
        actionButton(inputId="species",label = "species"),
        actionButton(inputId="genesG",label = "genes"),
        actionButton(inputId="genesCortrees",label = "Clusters genes"),
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
    trees=read.tree(input$trees$datapath, keep.multi = TRUE)
    sps = trees[[1]]$tip.label
    if (!is.null(sps)){
      vars <- all.vars(parse(text = sps))
    }
  })
  outliers <- reactive({
    input$Outlier
  })
  k <- reactive({
    input$k
  })
  distance <- reactive({
    input$choice
  })
  
  
###########-----------------------------------

  ##event = DETECTION OUTLIER ON
  
  observeEvent(outliers(), {
    if (outliers()== "on"){
    trees = read.tree(input$trees$datapath, keep.multi = TRUE)
      RES = Phylter(trees, distance = distance(), k=k(), thres = 0.5)
    output$outputPhylter1 <- renderText ({
      if(length(RES$Complete$outgn)!=0){
        gn=paste(RES$Complete$outgn, ";")
        append (gn, "outlier genes: ", after=0)
      }
      else{
        "no outlier gene detected"
      }
    })
    output$outputPhylter2 <- renderText ({  
      if(length(RES$Complete$outsp)!=0){
        sp=paste(RES$Complete$outsp, ";")
        append (sp, "outlier species: ", after=0)
      }
      else{
        "no outlier specie detected"
      }
    }) 
    output$outputPhylter3 <- renderText ({  
      if(!is.null(RES$CellByCell$outcell)){
        ce = paste(RES$CellByCell$outcell[,1],RES$CellByCell$outcell[,2], ";")
        append (ce, "outlier gène/species: ", after=0)
      }
      else{
        "no outlier cell detected"
      }
    })
    observeEvent(input$genesG, {
      trees = read.tree(input$trees$datapath, keep.multi = TRUE)
      dist=trees2dist(input$trees,distance())
      if(is.null(rownames(dist$res4Cmat$G))){
        if(!is.null(names(trees))){
          rownames(dist$res4Cmat$G)<-names(trees)
        }
        else{
          rownames(dist$res4Cmat$G)<-as.character(c(1:length(trees)))
        }
      }
      output$plot1 <- renderPlot({
        gn=paste(RES$Complete$outgn)
        G=dist$res4Cmat$G
        participant.colors <- rownames(G)
        if (length(RES$Complete$outgn)!=0){
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
          trees = read.tree(input$trees$datapath, keep.multi = TRUE)
          dist=trees2dist(input$trees,distance())
          if(is.null(rownames(dist$res4Cmat$G))){
            if(!is.null(names(trees))){
              rownames(dist$res4Cmat$G)<-names(trees)
            }
            else{
              rownames(dist$res4Cmat$G)<-as.character(c(1:length(trees)))
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
    dist = trees2dist(input$trees,distance())
    output$plot1 <- renderPlot({
    plot(hclust(as.dist(dist$res4Cmat$C)))
    })
  })
  observeEvent(input$species, {
    dist = trees2dist(input$trees,distance())
    FS=dist$res4Splus$F
    PF=dist$res4Splus$PartialF
    if(is.null(dimnames(PF)[[3]])){
      if(!is.null(names(input$trees))){
        dimnames(PF)[[3]] = names(input$trees)
      }
      else {
        dimnames(PF)[[3]] = as.character(c(1:length(input$trees)))
      }
    }
    observe({
      updateSelectInput(session, "selectSpecies", selected= "all", choices = append(outVar(),"all",after=0))
    })
    observeEvent(input$selectSpecies, {
        trees = read.tree(input$trees$datapath, keep.multi = TRUE)
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
       else{ 
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
      WR = trees2WR(input$trees,distance())
      plot.2WR(WR)
    })
  })
}

trees2dist <- function(trees, distance){
  trees = read.tree(trees$datapath, keep.multi = TRUE)
  if(is.null(names(trees))){names(trees)=as.character(c(1:length(trees)))}
  mat <- trees2matrices.Distatis(trees = trees, distance = distance)
  dist <- mat2Dist(mat)
  return(dist)
}

trees2WR <- function(trees, distance){
  dist <- trees2dist(trees = trees, distance = distance)
  WR=Dist2WR(dist)
  return(WR)
}

shinyApp(ui=ui,server =server)
