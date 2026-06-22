#### Style
shinyUI_taglist <- tagList(
				tags$style(HTML('* {font-family: "Calibri"}; ')),
				tags$style(HTML("input[type=\"number\"] {height: 30px;}")),
				tags$style(type = 'text/css', 
					'.navbar-default .navbar-brand {color: #ffffff; font-size:180%}',
					'.nav-item dropdown { font-size: 16px; font-size:140%; color: #ffffff}', 
					'.navbar .navbar-nav { font-size: 16px; font-size:140% }', 
					".nav-tabs {font-family: Calibri; font-size: 15px} ",
					".box {font-family: Calibri; font-size: 20px; word-wrap: break-word; } ",
					'.body {font-family: Calibri; font-size: 14px; word-wrap: break-word}',
					".navbar-default .navbar-brand {margin-left:-30px; margin-right:-30px; margin-top:-10px; padding-right:40px;}",
					".checkbox { margin-bottom: 0px; margin-left: 0px; margin-right: 0px; padding-left: 0px; padding-right: 0px;}",
					".table { margin-bottom: 0px; margin-left: 0px; margin-right: 0px; padding-left: 0px; padding-right: 0px;}"),
				tags$head(HTML('<link rel="icon", href="evergene.png", type="image/png" />')),
				tags$head(tags$style(HTML(".fileinput_2 {width: 0.1px;height: 0.1px;opacity: 0;overflow: hidden;position: absolute;z-index: -1;}")))
	)
				
