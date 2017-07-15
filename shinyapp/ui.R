library(shiny)

genes = as.character(read.table("../output/lof_enrichment_analysis/shiny_genes.tsv", 
                                header=FALSE)$V1)

shinyUI(pageWithSidebar(
  headerPanel('PTV Association Power'),
  sidebarPanel(
    selectInput("gene", "Gene", choices=genes, selectize=TRUE, selected="PCSK9"),
    numericInput("pD", "Disease prevalence", value=0.1, min=0, max=1, step=0.1),
    numericInput("RRAa", "Heterozygous relative risk", value=2),
    numericInput("nCase", "Number of cases", value=25000),
    numericInput("nControl", "Number of controls", value=25000),
    numericInput("alpha", "Type I error rate", value=2e-6),
    checkboxInput("unselected", label = "Unselected controls", value = FALSE)
  ),
  mainPanel(
    plotOutput('plot1')
  )
))