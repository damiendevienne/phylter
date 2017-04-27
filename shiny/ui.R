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
    tabPanel("visualize species on distatis compromise",
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
               The darker a couple gene/specie is, the further it is from it's reference position. You can zoom by brushing and doubleclicking on the plot"),
      plotOutput(outputId = "plot4", height = "1000px", dblclick="plot4dbclick", brush = brushOpts(id = "plot4_brush",resetOnNew = TRUE))
    ),
    tabPanel("visualize genes clusters",
      helpText("Clusterization of genes. Y axe is the correlation between genes"),
      plotOutput(outputId = "plot5", height = "1000px")
    ),
    tabPanel("visualize distances between species",
      helpText("Grey circles are the mean of distance between species for a particular species. Red circles are genes. 
               For a given gene, the distance between the specie of interest and an other one is given by the shape of the red line: 
               an acute angle mean an important distance between two species"),
      plotOutput(outputId = "plot6", height = "1000px")
    ),
    tabPanel("visualize distances between genes",
      helpText("For a given gene, each red line correspond to a species and each point on the red line correspond to 
               the distance between this specie and the other species? Grey circle is the mean distance between species for this given gene."),
      selectInput(inputId="Geneslist",label = "Select Genes",""),
      plotOutput(outputId = "plot7", height = "1000px")
    )
  )            
)
