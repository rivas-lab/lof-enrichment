library(shiny)
library(GeneticsDesign)

af.data = read.table("../output/lof_enrichment_analysis/shiny_data.tsv", 
                     sep="\t", header=TRUE, row.names=1)#, quote="\"")

shinyServer(function(input, output, session) {
  
  gene <- renderPrint(input$gene)
  pD <- renderPrint(input$pD)
  RRAa <- renderPrint(input$RRAa)
  nCase <- renderPrint(input$nCase)
  nControl <- renderPrint(input$nControl)
  alpha <- renderPrint(input$alpha)
  unselected <- renderPrint({ input$unselected })
  
  ratio = nControl / nCase
  
  power = sapply(d, function(x) GPC.default(pA=x, pD=pD, RRAa=RRAa, 
                                            RRAA=RRAa * 2, Dprime=1, 
                                            pB=pA, nCase=nCase, ratio=ratio, 
                                            alpha=alpha, unselected=unselected, 
                                            quiet=TRUE)$power)
  
  # Combine the selected variables into a new data frame
  selectedData <- reactive({
    iris[, c(input$xcol, input$ycol)]
  })

  clusters <- reactive({
    kmeans(selectedData(), input$clusters)
  })

  output$plot1 <- renderPlot({
    palette(c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3",
      "#FF7F00", "#FFFF33", "#A65628", "#F781BF", "#999999"))

    par(mar = c(5.1, 4.1, 0, 1))
    plot(selectedData(),
         col = clusters()$cluster,
         pch = 20, cex = 3)
    points(clusters()$centers, pch = 4, cex = 4, lwd = 4)
  })

})