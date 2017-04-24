ui<-fluidPage(
  theme = shinytheme("united"),
  wellPanel(
    fluidRow(
      column(6,
        fileInput(inputId = "trees", label = "add a tree file"),
        radioButtons(inputId="Outlier",label = "Detection", choiceNames = c("off","on"), choiceValues = c("off","on"), selected = "off")
      ),
      column(6,
        radioButtons(inputId="choice",label = "Method", choiceNames = c("nodal","patristic"), choiceValues = c("nodal","patristic"), selected = "patristic"),
        sliderInput(inputId="k", label = "select k", value = 3, min=1, max=6, round=FALSE, step=0.1)
      )
    ),
    fluidRow(
      textOutput(outputId= "outputPhylter1"),
      textOutput(outputId= "outputPhylter2"),
      textOutput(outputId= "outputPhylter3")
    )
  ),
  navbarPage("PhylteR",
    tabPanel("visualize some trees",
      column(3, offset = 0,
        textAreaInput(inputId="Treeslist",label = "Enter some genes"),
        bsPopover("Treeslist", "help", "Enter one or several genes separated by a comma", placement = "bottom", trigger = "hover", options = NULL),
        actionButton(inputId="SumitTrees",label="Visualize")
      ),
      mainPanel(
        plotOutput(outputId = "plot1", height = "1000px")
      )
    ),
    tabPanel("visualize genes by species",
      column(3, offset = 0,
        selectInput(inputId = "selectSpecies",label = "see a particular specie","all")
       ),
       mainPanel(
         helpText("Projection of the cross product of each gene matrice onto the compromise space of the distatis method"),
         plotOutput(outputId = "plot2", height = "1000px",dblclick = "BackClick2")
      )
    ),
    tabPanel("visualize genes",
      helpText("Plot of the between genes space. Genes with larger projections on the first axe are more similar to the other genes 
               than genes with smaller projections"),
      plotOutput(outputId = "plot3", height = "1000px")
    ),
    tabPanel("visualize 2WR matrices",
      helpText("The 2WR matrix is computed by calculating, for every species, the distance separating its position in each gene tree to its mean position. 
               The darker a couple gene/specie is, the further it is from it's reference position."),
      plotOutput(outputId = "plot4", height = "1000px")
    ),
    tabPanel("visualize genes clusters",
      plotOutput(outputId = "plot5", height = "1000px")
    ),
    tabPanel("visualize distances between species",
      plotOutput(outputId = "plot6", height = "1000px")
    ),
    tabPanel("visualize distances between genes",
      selectInput(inputId="Geneslist",label = "Select Genes",""),
      plotOutput(outputId = "plot7", height = "1000px")
    )
  )            
)
