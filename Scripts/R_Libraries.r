R_Libraries <- function(install_packs){
	if(install_packs == T){install.packages('XML')}else{suppressMessages(require(XML))}
	theurl <- "http://www.computerworld.com/article/2921176/business-intelligence/great-r-packages-for-data-import-wrangling-visualization.html"
	tables <- readHTMLTable(theurl)
	if(install_packs == T){
		install.packages(as.character(tables$cwsearchabletable$Package), repos = 'https://cran.mtu.edu/')
		devtools::install_github("muschellij2/matlabr")
		install.packages('R.matlab',repos = 'https://cran.mtu.edu/')
		install.packages("colorout")
	}
	# for(index in 1:length(as.character(tables$cwsearchabletable$Package))){suppressMessages(require(as.character(tables$cwsearchabletable$Package[index]), character.only = T))}
	suppressMessages(require(RPostgreSQL))
	suppressMessages(require(dplyr))
	suppressMessages(require(plyr))
	suppressMessages(require(readxl))
	suppressMessages(require(testthat))
	suppressMessages(require(dplyr))
	suppressMessages(require(csvread))
	suppressMessages(require(lubridate))
	suppressMessages(require(xlsx))
	suppressMessages(require(curl))
	suppressMessages(require(plyr))
	suppressMessages(require(XML))
	suppressMessages(require(RCurl))
	suppressMessages(require(data.table))
	suppressMessages(require(ggplot2))
	suppressMessages(require(TeachingDemos))
	suppressMessages(require(maps))
	# suppressMessages(require(mapdata))
	suppressMessages(require(maptools))
	suppressMessages(require(scales))
	suppressMessages(require(ggmap))
	suppressMessages(require(ggplot2))
	suppressMessages(require(grid))
	suppressMessages(require(GGally))
	suppressMessages(require(gridExtra))
	suppressMessages(require(jpeg))
	suppressMessages(require(R.matlab))
	suppressMessages(require(png))
	suppressMessages(require(shape))
	suppressMessages(require(DDHFm))
	suppressMessages(require(tiff))
	suppressMessages(require(dplyr))
	suppressMessages(require(igraph))
	# suppressMessages(require(devtool))
	suppressMessages(require(matlabr))
	suppressMessages(require(colorout))
	options(matlab.path = "/Applications/MATLAB_R2016b.app/bin")
}
