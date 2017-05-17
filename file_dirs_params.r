# Code from http://derekyves.github.io/2016/05/10/codeshare.html
# Load necessary libraries
setwd('MSP_Model/')
source('R_Libraries.r')
choice <- ifelse(readline(prompt = 'Install? y/n') == 'Y',TRUE,FALSE)
print(choice)
R_Libraries(choice) # After the first initial run this can be set to F
# Make a new environment:
fdirs <- new.env()

get_os <- function(){
	sysinf <- Sys.info()
	if (!is.null(sysinf)){
		os <- sysinf['sysname']
		if (os == 'Darwin')
			os <- "osx"
	} else {
		os <- .Platform$OS.type
		if (grepl("^darwin", R.version$os))
			os <- "osx"
		if (grepl("linux-gnu", R.version$os))
			os <- "linux"
	}
	tolower(os)
}
fdirs$computeros <- get_os()

## Declare the root project and data directory:
if(!grepl("windows", fdirs$computeros)){
  # Set R directories
  fdirs$home <- "~/MSP_Model/"
  fdirs$scrpdir <- "~/MSP_Model/Scripts/"
  fdirs$inpdatadir <- "~/MSP_Model/Input/Data/"
  fdirs$inpfigdir <- "~/MSP_Model/Input/Figures/"
  fdirs$outdatadir <- "~/MSP_Model/Output/Data/"
  fdirs$outfigdir <- "~/MSP_Model/Output/Figures/"

  # Set MATLAB Directories
  writeMat('~/MSP_Model/Scripts/Halibut/file_dir_params.mat',home_directory = fdirs$home, #cd(home_directory)
              script_dir = fdirs$scrpdir, #addpath(script_dir) % Script Directory
              input_figure_dir = fdirs$inpfigdir, #addpath(input_figure_dir) % Input Figure Directory
              input_data_dir = fdirs$inpdatadir, #addpath(input_data_dir) % Input Data Directory
              output_figure_dir = fdirs$outfigdir, #addpath(output_figure_dir) % Output Figure Directory
              output_data_dir = fdirs$outdatadir) #addpath(output_data_dir) % Output Data Directory)
}else{
  # Set R Directories
	root <- readline(prompt = "What is the name of the root drive? (Example: entering C would mean the model is located in drive C:/) ")
  fdirs$home <- paste0(root,"/MSP_Model/")
  fdirs$scrpdir <- paste0(root,"/MSP_Model/Scripts/")
  fdirs$inpdatadir <- paste0(root,"/MSP_Model/Input/Data/")
  fdirs$inpfigdir <- paste0(root,"/MSP_Model/Input/Figures/")
  fdirs$outdatadir <- paste0(root,"/MSP_Model/Output/Data/")
  fdirs$outfigdir <- paste0(root,"/MSP_Model/Output/Figures/")
  # Set MATLAB Directories
  writeMat('C:/MSP_Model/Scripts/Halibut/file_dir_params.mat',home_directory = fdirs$home, #cd(home_directory)
              script_dir = fdirs$scrpdir, #addpath(script_dir) % Script Directory
              input_figure_dir = fdirs$inpfigdir, #addpath(input_figure_dir) % Input Figure Directory
              input_data_dir = fdirs$inpdatadir, #addpath(input_data_dir) % Input Data Directory
              output_figure_dir = fdirs$outfigdir, #addpath(output_figure_dir) % Output Figure Directory
              output_data_dir = fdirs$outfigdir) #addpath(output_data_dir) % Output Data Directory)
}
# Load model
setwd(fdirs$home)
# Where is matlab?
input <- readline(prompt = "Where is MATLAB located? (Default location on OSX: /Applications/MATLAB_R2016b.app/bin/matlab ")
if(input == ""){
	print("Setting MATLAB location to /Applications/MATLAB_R2016b.app/bin/matlab")
	matlab_root <- "/Applications/MATLAB_R2016b.app/bin/matlab"
}else{
	print(paste0("Setting MATLAB location to ",input))
	matlab_root <- input
}
source(paste0(fdirs$scrpdir,'SCB_MSP_Model.r'))
