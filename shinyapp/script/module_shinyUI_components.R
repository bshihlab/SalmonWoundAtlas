#### module UI for making plotly plots (used across all plotly plots)
plotlyUI <- function(id, height ) {
	ns <- NS(id)
	tagList(
		plotlyOutput(ns("showPlot"), height = height)
	)
}
#### module UI for making base plots (used across all base plots)
plotlUI <- function(id, height ) {
	ns <- NS(id)
	tagList(
		plotOutput(ns("showPlot"), height = height)
	)
}
#### module UI for making base plots (used across all base plots)
buttonDownloadUI <- function(id, label) {
	ns <- NS(id)
	tagList(
		actionButton(
			inputId = ns("action"),
			label = label,
			icon = icon("download")
		)
	)
}


#### Panels
	
#### Main output panel 
### Subpanels within the main output panel
## Drop-down menus for selecting gene/cancer/pcx/pcy ----------
shinyUI_Panel_main1 <- fluidRow(
		conditionalPanel( condition = "output.panelStatus", 			
			column(width = 4, uiOutput("display_cancer_type")) %>% helper(type = "inline", title = "Row menu", content = help_rowMenu, buttonLabel="Close") ,
			column(width = 4, uiOutput("display_gene")),  
			column(width = 2, uiOutput("display_pcx")),  
			column(width = 2, uiOutput("display_pcy"))
		)
	)
		
## Gene contribution plot
shinyUI_Panel_main2 <- fluidRow(
		conditionalPanel( condition = "output.panelStatus", 
			#tags$h4("Gene contributions to principal components (PCs)"),
			plotlyUI("display_gene_contribution", height = "70px") %>% 
				helper(type = "inline", title = "Gene contribution to principal components (PCs)", content = help_gene_contribution, buttonLabel="Close")
		)
	)
## Principal component analysis-----------------------------------------
shinyUI_Panel_main3_tabPCA <- tabPanel("Principal component analysis (PCA)",
	fluidRow(
		conditionalPanel(condition = "output.panelStatus", 
			# Drop-down menu for sample annotation ----------
			column(width = 4, tags$h4("Main PCA Plots")) %>% 
			  helper(type = "markdown", title = "Additional Principal component analysis (PCA) plots", content = "pca1", buttonLabel="Close"),
			column(width = 4, uiOutput("display_pca_annotation")),
			column(width = 2, checkboxInput("pca_annotation_figLegend", "Hide figure legend", FALSE)),
			column(width = 2, checkboxInput("pca_annotation_na_rm", "Hide missing values", FALSE))
		)
	 ), # Fluid row

	fluidRow(
	   # Display PCA plot----------
	   splitLayout(#cellWidths =c("45%", "55%"),
			plotlyUI("display_pca_plot1", height = "400px"),
			plotlyUI("display_pca_plot2", height = "400px") 
		),
	), # Fluid row
	br(),
	fluidRow(
		 conditionalPanel(condition = "output.panelStatus", 
			helper(h4("Additional PCA plots"), type = "markdown", title = "Additional PCA plots", content = "pca2", buttonLabel="Close")
		)
	),# Fluid row
	## Additional PCA plots ---------------
	fluidRow(
		splitLayout(#cellWidths = c("45%", "55%"),
		# Display scree plot----------
			plotlyUI("display_scree_plot", height = "400px"),
			plotlyUI("display_geneAnnotation_plot", height = "400px")
		)
	), # Fluid row
	br(),

	##  Download buttons----------
	conditionalPanel(condition = "output.panelStatus", 
		fluidRow(
				helper(h4("Download"), type = "markdown", title = "Download", content = "pca3", buttonLabel="Close")
			
		),# Fluid row
		
		fluidRow(
			   column(width = 4, buttonDownloadUI("download_pdf_pca_button", "PDF (plots)")), 
			   column(width = 4, buttonDownloadUI("download_zip_pca_button", "ZIP (data and plots)"))
		), # Fluid row
	 ),
	 
	br(),
	fluidRow(
		 conditionalPanel(condition = "output.panelStatus", 
			helper(h4("Top gene contirbutions in each PC for the selected cancer type"), type = "inline", title = "Top gene contirbutions in each PC for the selected cancer type", content = help_gene_contribution_heatmap, buttonLabel="Close")
		)
	),# Fluid row
	
	fluidRow(
		conditionalPanel(condition = "output.panelStatus", 
		   # tile plot showing top contributing Genes
		   plotlyUI("display_tile_plot", height = "500px")
		)) # Fluid row
	) # Tab panel
## Survival analysis------------------------------
shinyUI_Panel_main3_tabSurvival <- tabPanel("Survival Analysis",
		# Survival input slider
		fluidRow(
			conditionalPanel( condition = "output.panelStatus", 
				# Drop-down menu for sample annotation ----------
				column(width = 1, tags$h4("")) %>% 
					helper(type = "markdown", title = "Survival analysis results", content = "survival1", buttonLabel="Close"),
				column(width = 2, numericInput("survival_percent_lower", "Low.exprs (%)", value = 33.3, min = 5, max = 95)) ,
				column(width = 2, numericInput("survival_percent_upper", "High.exprs (%)", value = 66.7, min = 5, max = 95)),
				column(width = 3, uiOutput("display_survival_event")),
				column(width = 4, uiOutput("display_survival_annotation")),
			),
		), # Fluid row
		fluidRow(
			conditionalPanel( condition = "output.panelStatus", 
				# Drop-down menu for sample annotation ----------
				column(width = 6, tags$p("")) ,
				column(width = 3, checkboxInput("survival_annotation_figLegend", "Hide figure legend", FALSE)),
				column(width = 3, checkboxInput("survival_annotation_na_rm", "Hide missing values", FALSE))
			)														
		),
		 # Survival main
		fluidRow(
			 splitLayout(
					plotlUI("display_survival_plot1", height = "400px"),
					plotlyUI("display_survival_plot2", height = "400px"))
		), # Fluid row
		br(),
		 # Download button (reactive)
		conditionalPanel(condition = "output.panelStatus", 
			fluidRow(
					helper(h4("Download"), type = "markdown", title = "Download", content = "survival2", buttonLabel="Close")
				
			),# Fluid row
		fluidRow(
			column(width = 5, buttonDownloadUI("download_pdf_survival_button", "PDF (plots)")), 
			column(width = 5, buttonDownloadUI("download_zip_survival_button", "Zip (data and plots)"))
		#column(width = 5, uiOutput("download_csv_survival_button"))
		) # Fluid row
		)
	) # Tab panel
	
## Correlation analysis------------------------------
shinyUI_Panel_main3_tabCorr <- tabPanel("Correlation",
	# Survival input slider
	fluidRow(
		conditionalPanel( condition = "output.panelStatus", 
			# Drop-down menu for sample annotation ----------
			#column(width = 1, tags$h4("")) 
			column(width = 5, uiOutput("display_correlation_annotation1"))%>% 
				helper(type = "markdown", title = "Correlation", content = "correlation", buttonLabel="Close"),
			column(width = 5, uiOutput("display_correlation_annotation2")),
			#column(width = 5, pickerInput("display_correlation_annotation2", "Y-axis (up to 9)", options = list('actions-box' = TRUE), list_of_cancer_types, multiple = TRUE)),
			#pickerInput("cancer_type_list", "Cancer Types", options = list('actions-box' = TRUE), list_of_cancer_types, multiple = TRUE) 
			column(width = 2, checkboxInput("correlation_log_gene_exp", "log gene expression", TRUE)))
			), # Fluid row
	 # Survival main
	 fluidRow(
		splitLayout(cellWidths = c("80%", "20%"),
			plotlyUI("display_corr_plot", height = "500px"),
			plotlyOutput("display_corr_plot0", height = "500px"))
			),# Fluid row
		fluidRow(
			column(width = 5, buttonDownloadUI("download_pdf_correlation_button", "PDF (plots; selected x-axis vs all)")), 
			column(width = 5, buttonDownloadUI("download_zip_correlation_button", "Zip (data and plots; selected x-axis vs all)"))
		) # Fluid row
	) # Tab panel
	
	
## Side panel for cancer/gene input
shinyUI_Panel_inupt <- sidebarPanel( 
		width = 3,
		# User input
		prettySwitch(inputId = "user_input_show", label = "User input", fill = TRUE, status = "primary", value=FALSE) %>% 
			helper(type = "markdown", title = "User input data", content = "userinput", buttonLabel="Close"),

		# Cancer type input			
		conditionalPanel( condition = "!input.user_input_show", 
			pickerInput("cancer_type_list", "Cancer Types", options = list('actions-box' = TRUE), list_of_cancer_types, multiple = TRUE)  %>% 
				helper(type = "inline", title = "Cancer types", content = help_cancer_input, buttonLabel="Close")
		),

		conditionalPanel( condition = "input.user_input_show", 
			fileInput("upload_count", "Data", multiple = FALSE, accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv")),
			fileInput("upload_annotation", "Sample annotation", multiple = FALSE, accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv")),
			fileInput("upload_gene", "Gene annotation (optional)", multiple = FALSE, accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv")),
			fileInput("upload_survival", "Survival data (optional)", multiple = FALSE, accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv")),
			#fileInput_alt("upload_count", " Data", progress=TRUE),
			#fileInput_alt("upload_annotation", " Sample annotation ", progress=TRUE),
			#fileInput_alt("upload_gene", " Gene annotation ", progress=TRUE),
			#fileInput_alt("upload_survival", " Survival data ", progress=TRUE),
			checkboxInput("upload_log_transform", "Add log transformation", FALSE),
			checkboxInput("upload_raw_count", "Raw count", FALSE),
		),
		# Tag
		#tags$h4("OR"),
		#tags$p("Input data (PCA and correlation only)"),
		# User input
		#fileInput("upload_count_mx", label=NULL),
		#fileInput("upload_sample_annotation", label=NULL),

		#textOutput("computed"),
		# Gene input
		textAreaInput("gene_string", "Genes", rows = 5) %>% 
			helper(type = "inline", title = "Genes", content = help_gene_input, buttonLabel="Close"),
		#tags$p("Maximum number of gene-cancer combination is 100."),
		
		
		# Reactive error message
		textOutput("invalid_gene_message_display"),
		textOutput("input_project_update"),
		
		
		# Compute button
		analyseButton("button_compute"),
		br(),
		exampleButton("example_case1", label= "Example 1") %>% 
		  helper(type = "inline", title = "Example 1: esophageal carcinoma", content = help_example1, buttonLabel="Close"), 
		exampleButton("example_case2", label= "Example 2") %>% 
		  helper(type = "inline", title = "Example 2: breast invasive carcinoma", content = help_example2, buttonLabel="Close"), 
		br(),
										
		br(),
		tags$p("Please interpret all results in evergene with caution as they can be confounded by many variables that are not controlled for in these analysis. Consider consulting a statistician."),
		br()

	) # Sidebar panel 
