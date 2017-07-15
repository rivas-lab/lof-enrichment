library(shiny)
library(GeneticsDesign)
library(ggplot2)

af.data <- read.table("../output/lof_enrichment_analysis/shiny_data.tsv", 
                     sep="\t", header=TRUE, row.names=1)#, quote="\"")
theme_set(theme_bw(base_size=18))
a = colnames(af.data)[seq(1, dim(af.data)[2], 3)]
xlabels <- as.factor(sapply(a, function(x) unlist(strsplit(x, split="_"))[3]))

shinyServer(function(input, output, session) {
  
  afs <- reactive({af.data[input$gene, ]})
  ratio <- reactive({input$nControl / input$nCase})
  
  power.vals <- reactive({sapply(afs(), function(x) GPC.default(
    pA=x, pD=input$pD, RRAa=input$RRAa, RRAA=input$RRAa * 2, 
    Dprime=1, pB=x, nCase=input$nCase, ratio=ratio(), alpha=input$alpha, 
    unselected=input$unselected, quiet=TRUE)$power)
  })

  power.df <- reactive({data.frame(
    x = seq(1, length(power.vals()) / 3, 1),
    y = power.vals()[seq(1, length(power.vals()), 3)],
    ymin = power.vals()[seq(2, length(power.vals()), 3)],
    ymax = power.vals()[seq(3, length(power.vals()), 3)],
    xlabels = xlabels
  )})

  output$plot1 <- renderPlot({
    ggplot(power.df(), aes(x = x,y = y)) + 
      geom_point(size=4) + 
      geom_errorbar(aes(ymin = ymin, ymax = ymax), width=0.25) + 
      xlab("Populations") + 
      ylab("Power") + 
      scale_x_continuous(breaks=seq(7), labels = power.df()$xlabels)
  })

})