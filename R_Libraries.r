Useful_R_packages <- function(install_packs = T, print_names = T){
	suppressMessages(require(XML))
	theurl <- "http://www.computerworld.com/article/2921176/business-intelligence/great-r-packages-for-data-import-wrangling-visualization.html"
	tables <- readHTMLTable(theurl)
	if(install_packs){
		install.packages(as.character(tables$cwsearchabletable$Package), repos = 'https://cran.mtu.edu/')
		devtools::install_github("hadley/readxl")
	}
	for(index in 1:length(as.character(tables$cwsearchabletable$Package))){
		if(print_names){
			print(as.character(tables$cwsearchabletable$Package[index]))
		}
		suppressMessages(require(as.character(tables$cwsearchabletable$Package[index]), character.only = T))
	}
	suppressMessages(require(RPostgreSQL))
	suppressMessages(require(dplyr))
	suppressMessages(require(plyr))
	suppressMessages(require(readxl))
	suppressMessages(require(testthat))
}
