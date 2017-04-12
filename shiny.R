library(shiny)
library(shinythemes)
library(shinyHeatmaply)
source("/home/aurore/Documents/Phylter/pmcoa.R")
source("/home/aurore/Documents/Phylter/PhylteR.R")

#trees = read.tree("/home/aurore/Documents/Phylter/viz/Aguileta-et-al-2008_TREES.txt", keep.multi = TRUE)
#trees=read.tree(file="/home/aurore/Documents/Phylter/trees/rose/test.phy", keep.multi = TRUE)
ui<-fluidPage(
  theme = shinytheme("flatly"),
  fluidRow(
    column(3,
      wellPanel(
        titlePanel("Data input"),
        fileInput(inputId = "trees", label = "add a tree file")
      ),
      wellPanel(
        titlePanel("Detection"),
        checkboxGroupInput(inputId="choice",label = "Method", choiceNames = c("nodal","patristic"), choiceValues = c("nodal","patristic")),
        actionButton(inputId="Outlier",label = "Outlier Detection"),
        sliderInput(inputId="k", label = "select k", value = 2, min=1, max=6,round=FALSE,step=0.1)
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
    column(9,
      plotOutput(outputId = "plot1", height = "1000px")   
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
  
  observeEvent(input$WR, {
    output$plot1 <- renderPlot({
    #dist = trees2dist(input$trees,input$choice)
    WR = trees2WR(input$trees,input$choice)
    plot.2WR(WR)
    })
  })
  observeEvent(input$genesCortrees, {
    dist = trees2dist(input$trees,input$choice)
    output$plot1 <- renderPlot({
    plot(hclust(as.dist(dist$res4Cmat$C)))
    })
  })
  observeEvent(input$species, {
    dist = trees2dist(input$trees,input$choice)
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
      updateSelectInput(session, "selectSpecies", selected= "all", choices = append(all.vars(parse(text = sps)),"all"))
    })
    observeEvent(input$selectSpecies, {
        trees = read.tree(input$trees$datapath, keep.multi = TRUE)
        S=input$selectSpecies
        print(S)
       # if (S == "all"){
        #  output$plot1 <- renderPlot({
         #   par(mfrow=c(ceiling(length(trees[[1]]$tip.label)/4),4))
         #   par(mar=c(2,2,2,2))
         #   part.design <- diag(dim(PF)[3])
          #  participant.colors <- as.matrix(createColorVectorsByDesign(part.design)$oc)
         #   for (sp in 1:length(trees[[1]]$tip.label)){
          #    to.plot <- t(PF[sp, , ])
           #   center.point <- FS[sp, c(1, 2)]
          #    center.rep <- matrix(center.point, dim(PF)[3], 2, byrow = TRUE)
           #   bound.mat <- rbind(center.rep, to.plot[, c(1, 2)])
           #   bound.mat <- bound.mat[as.vector(t(matrix(seq(1, nrow(bound.mat)), ncol = 2))), ]
            #  plot(to.plot, main = dimnames(PF)[[1]][sp], cex.main= 1.5, col = participant.colors, cex=3.2)
            #  points(bound.mat, type = "l", lty = 1, lwd = 1, col = "grey70")
            #  text(to.plot, labels = rownames(t(PF[sp, , ])))
            #}
         # })
        #}
       #else{ 
       #  output$plot1 <- renderPlot({
       #   part.design <- diag(dim(PF)[3])
        #  participant.colors <- as.matrix(createColorVectorsByDesign(part.design)$oc)
       #   to.plot <- t(PF[S, , ])
        #  center.point <- FS[S, c(1, 2)]
        #  center.rep <- matrix(center.point, dim(PF)[3], 2, byrow = TRUE)
        #  bound.mat <- rbind(center.rep, to.plot[, c(1, 2)])
        #  bound.mat <- bound.mat[as.vector(t(matrix(seq(1, nrow(bound.mat)), ncol = 2))), ]
        #  plot(to.plot, main = dimnames(PF)[[1]][S], cex.main= 1.5, col = participant.colors, cex=3.2)
        #  points(bound.mat, type = "l", lty = 1, lwd = 1, col = "grey70")
        #  text(to.plot, labels = rownames(t(PF[S, , ])))
      })
      
    })
  observeEvent(input$genesG, {
    trees = read.tree(input$trees$datapath, keep.multi = TRUE)
    dist=trees2dist(input$trees,input$choice)
   if(is.null(rownames(dist$res4Cmat$G))){
      if(!is.null(names(trees))){
        rownames(dist$res4Cmat$G)<-names(trees)
      }
     else{
        rownames(dist$res4Cmat$G)<-as.character(c(1:length(trees)))
     }
  }
    output$plot1 <- renderPlot({
      participant.colors <- color.scale(dist$res4Cmat$G[,1],cs1=c(1,1,0),0)
      plot(dist$res4Cmat$G, col = participant.colors, cex = 3.5)
      text(dist$res4Cmat$G, labels = rownames(dist$res4Cmat$G))
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
