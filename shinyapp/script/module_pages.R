### each of the main pages
# Analysis
pages_analysis <- tabPanel("Analysis", 
						page_sidebar(
							# Sidebar panel for inputs ----
							sidebar = sidebar(
							# Input: Slider for the number of bins ----
							sliderInput(
							  inputId = "point_size",
							  label = "Plot point size",
							  min = 0.1,
							  max = 3,
							  value = 0.5
							),
							# genes input
							textAreaInput(
							  inputId = "input_genes",
							  label = "Input genes"
							) %>% 
							  helper(
							        "question-circle",
							        type = "inline", 
							         title = "Gene input", 
							         content ="Please enter up to 30 genes (1 gene per line). This can be Ensembl gene ID for Atlantic salmon or gene names."),
							# button
							actionButton("submit_button", "Submit"),
							textOutput("msgTxt")
							),
							
							# Input: select gene menu
							selectInput(
							  inputId = "display_gene",
							  label = "Gene",
							  choices = c("Please enter input genes on the left.")
							),
							
							input_switch("splitTreat", "Split by treatment"), 

							# Output: plot ----
							plotlyOutput("dotplot",height = "400px", width = "550px")
							#plotlyUI("display_pca_plot1", height = "400px"),
							)
						)
# About
pages_about <- tabPanel("About",

	) # Tab panel

# Contact us
pages_contactUs <- tabPanel("Contact us")

# Source page
pages_source <- tabPanel("Source code", 
	HTML('<script src="https://emgithub.com/embed-v2.js?target=https%3A%2F%2Fgithub.com%2Fellaintheclouds%2FNWCR%2Fblob%2Fmain%2FPCA_function.R&style=default&type=code&showBorder=on&showLineNumbers=on&showFileMeta=on&showFullPath=on&showCopy=on"></script>'))
	
# Documentation

#pages_documentation <- tabPanel("Home", href = "https://google.com", style = "color:white;")
pages_documentation <- tabPanel("Documentation", tags$iframe(style="height:800px; width:100%; scrolling=yes", src="documentation.pdf"))
	