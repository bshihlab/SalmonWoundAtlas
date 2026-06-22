### Combining the composit tabs/panels within in the main output panel
## Main3 (tabs within main output panel)
shinyUI_Panel_main3 <- tabsetPanel(type = "tabs", 
						shinyUI_Panel_main3_tabPCA, 
						shinyUI_Panel_main3_tabSurvival, 					
						shinyUI_Panel_main3_tabCorr) 

## Everything together
shinyUI_Panel_main <- mainPanel(			
					shinyUI_Panel_main1, 
					shinyUI_Panel_main2,
					shinyUI_Panel_main3			
				  ) # Main panel