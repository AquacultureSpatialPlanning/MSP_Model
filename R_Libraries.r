R_Libraries <- function(install_packs){
	if(install_packs == T){install.packages('XML')}else{suppressMessages(require(XML))}
	# theurl <- "http://www.computerworld.com/article/2921176/business-intelligence/great-r-packages-for-data-import-wrangling-visualization.html"
	# tables <- readHTMLTable(theurl)
	packs <- c('RPostgreSQL','dplyr','tidyr','plyr','readxl','testthat',
						'csvread','lubridate','xlsx','curl','XML','RCurl',
						'data.table','ggplot2','TeachingDemos','maps','maptools',
						'scales','ggmap','GGally','gridExtra','jpeg','R.matlab','png',
						'shape','DDHFm','tiff','igraph','devtools','matlabr','colorout',
						'XLConnect','knitr','gtools','TeachingDemos','grDevices','svg')
	if(install_packs == T){
		# install.packages(as.character(tables$cwsearchabletable$Package), repos = 'https://cran.mtu.edu/')
		devtools::install_github("muschellij2/matlabr")
		install.packages('R.matlab',repos = 'https://cran.mtu.edu/')
		install.packages("colorout")
		install.packages(packs,repos = 'https://cran.mtu.edu/')
	}
	lapply(packs,require,character.only = TRUE)
	library(grid)
	options(matlab.path = "/Applications/MATLAB_R2016b.app/bin")
}
