# Load necessary libraries
this.dir <- dirname(sys.frame(1)$ofile)
# Set R directories
home <- file.path(paste0(this.dir))
scrpdir <- file.path(paste0(this.dir,'/Scripts/'))
inpdatadir <- file.path(paste0(this.dir,'/Input/Data/'))
inpfigdir <- file.path(paste0(this.dir,'/Input/Figures/'))
outdatadir <- file.path(paste0(this.dir,'/Output/Data/'))
outfigdir <- file.path(paste0(this.dir,'/Output/Figures/'))

source(file.path(paste0(home,'/R_Libraries.r')))
choice <- ifelse(readline(prompt = 'Install? y/n ') == 'Y',TRUE,FALSE)
print(choice)
R_Libraries(choice)
# Set MATLAB Directories
writeMat(file.path(paste0(scrpdir,'/Halibut/file_dir_params.mat')),home_directory = home, #cd(home_directory)
						script_dir = scrpdir, #addpath(script_dir) % Script Directory
						input_figure_dir = inpfigdir, #addpath(input_figure_dir) % Input Figure Directory
						input_data_dir = inpdatadir, #addpath(input_data_dir) % Input Data Directory
						output_figure_dir = outfigdir, #addpath(output_figure_dir) % Output Figure Directory
						output_data_dir = outdatadir) #addpath(output_data_dir) % Output Data Directory)
setwd(home)
# Where is matlab?
input <- readline(prompt = "Where is MATLAB located? (Default location on OSX: /Applications/MATLAB_R2016b.app/bin/matlab ")
if(input == ""){
	print("Setting MATLAB location to /Applications/MATLAB_R2016b.app/bin/matlab")
	matlab_root <- "/Applications/MATLAB_R2016b.app/bin/matlab"
}else{
	print(paste0("Setting MATLAB location to ",input))
	matlab_root <- input
}
source(paste0(scrpdir,'SCB_MSP_Model.r'))
