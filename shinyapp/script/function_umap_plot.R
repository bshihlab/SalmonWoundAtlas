## umap plot
umapPlot <- function(id, display_cancer_type, display_gene, display_pcx, display_pcy, pca_annotation_na_rm, rv) {
  moduleServer( id,
                function(input, output, session) {
                  output$showPlot <-   output$dotplot <- renderPlotly({
                    
                    validate(need(display_cancer_type(), message = FALSE), 
                             need(display_gene(), message =  FALSE))
                    
                    # get current plot data
                    plot_df <- data.frame(exprs = rv$current_exprs_df[,input$display_gene])
                    #plot_df <- data.frame(exprs = exprs_df[,input$display_gene])
                    plot_df <- cbind(plot_df, cell_embedding, cell_metadata)
                    plot_df <- plot_df[order(plot_df$exprs),]
                    colnames(plot_df)[1] <- "exprs"
                    
                    # q95 max cut off
                    q95_cutoff <- quantile(plot_df$exprs[plot_df$exprs!=0],0.95)
                    plot_df$exprs <- ifelse(plot_df$exprs > q95_cutoff, q95_cutoff, plot_df$exprs)
                    
                    
                    #x    <- faithful$waiting
                    # make plotly plot
                    p <- plot_ly(data = plot_df, 
                                 type = "scatter",
                                 mode = "markers",
                                 x = ~UMAP_1, 
                                 y = ~UMAP_2,
                                 color = ~exprs, 
                                 text = ~named_celltypes,
                                 marker = list(size = input$point_size,
                                               reversescale=T))
                    ggplotly(p) 
                    
                  })
                  
                }
  )
}
