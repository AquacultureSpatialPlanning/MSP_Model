Requirements
1.) This folder must be put into the root folder (i.e. '~' in OSX and C:/ in Windows)
2.) Both MATLAB and R must be installed
3.) MATLAB must be installed with the following Toolbox Mapping Toolbox,Bioinformatics,Parallel Optimization

Instructions
1.) Clone repository into root directory
2.) Open terminal(osx) or command prompt (windows)

OSX Instructions
1.) If homebrew is not currently on the machine follow these instructions to do so https://www.moncefbelyamani.com/how-to-install-xcode-homebrew-git-rvm-ruby-on-mac/
2.) Once homebrew is installed follow these instuctions to download R https://rud.is/b/2015/10/22/installing-r-on-os-x-100-homebrew-edition/
3.) In terminal type, c ~, to go to the root directory
4.) Launch R in terminal
5.) Type source('~/MSP_Model/file_dir_params.r') and then follow the prompts

Windows Instrucions
1.) Follow these instructions to download R for windows https://cran.r-project.org/bin/windows/base/
2.) Launch R from the command line or in RStudio
3.) Run the command source('C:/MSP_Model/file_dir_params.r') (C:/ could be any drive as long as the folder MSP_Model is in the root directory)

Note: 
To run the code that produces the filter and seed plans:
-Open Matlab
-Open and run the file 'Filter_Seed_v6_BrayCurtis_publish_v1.m'
