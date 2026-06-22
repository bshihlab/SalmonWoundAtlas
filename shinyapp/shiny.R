# source functions and modules
source("script/load_libraries.R")
source("script/module_pages.R")


## read in data
# exprs_data <- readMM("shinyapp/data/h5ad_data.mx")
# cell_embedding <- read.csv("shinyapp/data/h5ad_embedding.csv")
# cell_metadata <- read.csv("shinyapp/data/h5ad_annotation_brief.csv")
# gene_df <- read.csv("shinyapp/data/gene.csv")
# gene_orthlog <- read.csv("shinyapp/data/orth_zebsalhum.csv")

#exprs_data <- readMM("../analysis/salmon_wound_atlas/rds/h5ad_data.mx")
cell_embedding <- read.csv("data/h5ad_embedding.csv")
cell_metadata <- read.csv("data/h5ad_annotation_brief.csv")
gene_df <- read.csv("data/gene.csv")

colnames(cell_embedding) <- c("UMAP_1", "UMAP_2")



#exprs_df <- read.csv("shinyapp/data/test_df.csv")
#exprs_df <- read.csv("data/test_df.csv")
#exprs_df <- reactiveValues(exprs_df = exprs_df)

# q95 max cut off
#q95_cutoff <- quantile(plot_df$exprs[plot_df$exprs!=0],0.95)
#plot_df$exprs <- ifelse(plot_df$exprs > q95_cutoff, q95_cutoff, plot_df$exprs)


## Define UI for app that draws a histogram ----
# ui <- page_sidebar(
#   # App title ----
#   title = "Hello Shiny!",
#   # Sidebar panel for inputs ----
#   sidebar = sidebar(
#     # Input: Slider for the number of bins ----
#     sliderInput(
#       inputId = "point_size",
#       label = "Plot point size",
#       min = 1,
#       max = 10,
#       value = 2
#     ),
#     # genes input
#     textAreaInput(
#       inputId = "genes",
#       label = "Genes"
#     ),
#     # button
#     actionButton("submit_button", "Submit")
#   ),
# 
# 
#   # Input: select gene menu
#   selectInput(
#     inputId = "display_gene",
#     label = "Gene",
#     choices = c("ENSSSAG00000042105", 
#                 "ENSSSAG00000046009", 
#                 "ENSSSAG00000069316", 
#                 "ENSSSAG00000050783")
#     ),
#     
#     # Output: Histogram ----
#     plotlyOutput("dotplot",height = "500px")
# )
ui <- fluidPage(theme = shinytheme("sandstone"),
                #shinyUI_taglist,
                navbarPage(
                  title="Salmon wound atlas",
                  pages_analysis,
                  #pages_about
                  #pages_documentation,
                  #pages_contactUs,
                  #pages_source
                )
)

## Define server logic required to draw a histogram ----
server <- function(input, output, session) {
  observe_helpers()
  rv <- reactiveValues(
                      fst_idx = NULL,
                       exprs_data = NULL,
                       current_exprs_df = NULL,
                       matched_genelist = NULL  ,
                       display_message = "")
  observeEvent(input$submit_button, {
    ## load data if not already loaded
    showPageSpinner(caption = "Loading, please wait...")
     if(is.null(rv[["fst_idx"]])){
       rv[["fst_idx"]] <- read_fst("data/exprs_rownum_indicate.fst")
     }
    ## Find the corresponding genes
    formatted_gene_list <- gsub(" ", "", input$input_genes)
    formatted_gene_list <- gsub(",", "\n", input$input_genes)
    formatted_gene_list <- strsplit(formatted_gene_list, "\n")[[1]]

    formatted_gene_list <- toupper(formatted_gene_list)
    formatted_gene_list <- formatted_gene_list[formatted_gene_list != ""]
  
    # cap input genes to 30 genes
    max_gene <- ifelse(length(formatted_gene_list) > 30, 30, length(formatted_gene_list))
    formatted_gene_list <- formatted_gene_list[1:max_gene]
    
    # Get matching salmon gene ids and their display gene names
    matched_sid <-  gene_df$salmon.geneid[gene_df$salmon.geneid %in% formatted_gene_list]
    matched_s <-  gene_df$salmon.geneid[toupper(gene_df$salmon.genename) %in% formatted_gene_list]
    matched_z <-  gene_df$salmon.geneid[toupper(gene_df$zebrafish.genename) %in% formatted_gene_list]
    matched_h <-  gene_df$salmon.geneid[toupper(gene_df$human.genename) %in% formatted_gene_list]
    matched_geneid <- unique(c(matched_sid, matched_s, matched_z, matched_h))
    gene_df_matched <- gene_df[gene_df$salmon.geneid %in% matched_geneid,]

    if(nrow(gene_df_matched)>0){
      gene_df_matched$combined_names <- paste0(gene_df_matched$gene_name_base, "_", gene_df_matched$salmon.geneid)
      gene_df_matched <- gene_df_matched[order(gene_df_matched$combined_names),]
      
      # keep a maximum of 30 genes
      #max_gene <- ifelse(nrow(gene_df_matched)>30, 30, nrow(gene_df_matched))
      #gene_df_matched <- gene_df_matched[1:max_gene,]
      rv[["matched_genelist"]] <- gene_df_matched$combined_names
      
      # get the row ids for the matching gene ids
      matched_row_num <- match( gene_df_matched$salmon.geneid, gene_df$salmon.geneid)
      print(gene_df_matched)
      current_exprs_list <- list()
      current_geneid_list <- vector()
      gene_count=1
      ###### need to loop through the 2 fst files, 
      # but not done this here since sharing through SCP upload 
      # would be better than making a shinyapp work
      for(fst_idx_row in matched_row_num){
        rowRange_df <- rv[["fst_idx"]][fst_idx_row,]
        data_df <- read_fst("data/exprs.fst", 
                            from = rowRange_df$rowMin, 
                            to = rowRange_df$rowMax)
        data_df$i <- 1
        data_df <- sparseMatrix( i = data_df$i, j = data_df$j, x = data_df$x, dims = c(1, 40127))
        current_exprs_list[[gene_count]] <- as.data.frame(t(as.matrix(data_df)))
        gene_count = gene_count+1
      }
      current_exprs_df <- do.call(cbind, current_exprs_list)
      colnames(current_exprs_df) <- gene_df_matched$combined_names
      rv[["current_exprs_df"]] <- current_exprs_df
      updateSelectInput(inputId = "display_gene", choices = rv[["matched_genelist"]])
      rv[["display_message"]] <- ""
      
    } else {
      rv[["display_message"]] <- "Gene input not found."
    }
    hidePageSpinner()
  })

  #### render text
  output$msgTxt <- renderText(rv[["display_message"]])
  
  #### Plot  
  output$dotplot <- renderPlotly({
    
    validate(need(rv$current_exprs_df, message = FALSE), 
             need(input$display_gene, message =  FALSE),
             need(input$display_gene %in% colnames(rv$current_exprs_df), message=FALSE))
    
    # get current plot data
    current_gene <- input$display_gene
    plot_df <- data.frame(exprs = rv$current_exprs_df[,current_gene])
    #plot_df <- data.frame(exprs = exprs_df[,input$display_gene])
    plot_df <- cbind(plot_df, cell_embedding, cell_metadata)
    plot_df <- plot_df[order(plot_df$exprs),]
    colnames(plot_df)[1] <- "exprs"
    
    # q95 max cut off
    q95_cutoff <- quantile(plot_df$exprs[plot_df$exprs!=0],0.95)
    plot_df$exprs <- ifelse(plot_df$exprs > q95_cutoff, q95_cutoff, plot_df$exprs)
    
    
    #x    <- faithful$waiting
    # make plotly plot
    p <- ggplot(data= plot_df,
                aes(x=UMAP_1, y=UMAP_2, colour=exprs, text=named_celltypes)) + 
          geom_point(size = input$point_size) +
          theme_bw() +
          scale_colour_viridis(option = "A",direction = -1)
    
    if(input$splitTreat){
      p <- p + facet_wrap(~treatment_group)
    }
    
    ggplotly(p, tooltip = c("named_celltypes"))
    
  })
  
}

shinyApp(ui = ui, server = server)

