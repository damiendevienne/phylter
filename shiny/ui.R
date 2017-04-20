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
        textAreaInput(inputId="Treeslist",label = "Enter some genes"),
        actionButton(inputId="SumitTrees",label="Visualize")
      ),
      wellPanel(
        actionButton(inputId="species",label = "SpeciesByGenes"),
        selectInput(inputId = "selectSpecies",label = "see a particular specie","")
      ),
      wellPanel(
        actionButton(inputId="species2",label = "Species"),
        actionButton(inputId="genesG",label = "Genes"),
        textAreaInput(inputId="Geneslist",label = "Enter some genes"),
        actionButton(inputId="genes2",label = "Genes2"),
        actionButton(inputId="genesCortrees",label = "Genes Clusters"),
        actionButton(inputId="WR",label = "WR")
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

