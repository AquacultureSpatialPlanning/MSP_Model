R_Libraries <- function(install_packs){
<<<<<<< HEAD
	suppressMessages(require(XML))
=======
	if(install_packs == T){install.packages('XML')}else{suppressMessages(require(XML))}
>>>>>>> e9a6acd68e4a0fd08b18829858208289c9350661
	theurl <- "http://www.computerworld.com/article/2921176/business-intelligence/great-r-packages-for-data-import-wrangling-visualization.html"
	tables <- readHTMLTable(theurl)
	if(install_packs == T){
		install.packages(as.character(tables$cwsearchabletable$Package), repos = 'https://cran.mtu.edu/')
<<<<<<< HEAD
		devtools::install_github("hadley/readxl")
		devtools::install_github("muschellij2/matlabr")
=======
		devtools::install_github("muschellij2/matlabr")
		install.packages('R.matlab',repos = 'https://cran.mtu.edu/')
>>>>>>> e9a6acd68e4a0fd08b18829858208289c9350661
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
<<<<<<< HEAD
	suppressMessages(require(mapdata))
=======
	# suppressMessages(require(mapdata))
>>>>>>> e9a6acd68e4a0fd08b18829858208289c9350661
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
<<<<<<< HEAD
	suppressMessages(require(devtool))
	suppressMessages(require(matlabr))
	suppressMessages(require(colorout))
	suppressMessages(require(matlabr))
=======
	# suppressMessages(require(devtool))
	suppressMessages(require(matlabr))
	suppressMessages(require(colorout))
>>>>>>> e9a6acd68e4a0fd08b18829858208289c9350661
	options(matlab.path = "/Applications/MATLAB_R2016b.app/bin")
}
