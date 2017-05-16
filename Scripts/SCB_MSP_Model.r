<<<<<<< HEAD
# NOTE: MATLAB requires 'Mapping Toolbox', 'Bioinformatics', 'Parallel Optimization'
=======
# MATLAB requires 'Mapping Toolbox', 'Bioinformatics', 'Parallel Optimization'
>>>>>>> e9a6acd68e4a0fd08b18829858208289c9350661
# Set current working directory as a string
wkdir <- getwd()
# Load necessary R libraries. For the function R_Libraries, enter T if this is the first time running the model. This will
# install all of the necessary libraries and load them into the
# current workspace.
source(paste0(wkdir,'/MSP_Model/Scripts','/R_Libraries.r'))
<<<<<<< HEAD
R_Libraries(F) # After the first initial run this can be set to F
=======
R_Libraries(FALSE) # After the first initial run this can be set to F
# Install R markdown
install.packages("knitr",repos = 'https://cran.mtu.edu/')
library(knitr)
>>>>>>> e9a6acd68e4a0fd08b18829858208289c9350661
# Set global variables
n.sector <- 7 # Number of sectors
epsilon <- 0.2 # Stepsize of sector weights
t <- 10 # Time Horizon
r <- 0.05 # Discount rate
<<<<<<< HEAD
paste0(
=======

>>>>>>> e9a6acd68e4a0fd08b18829858208289c9350661
# Read sector data
sector_data.df <- read.csv(paste0(wkdir,'/MSP_Model/Input/Data/SeaGrant_data_complete_2015.csv'))
fulldomain <- sector_data.df$TARGET_FID # Model domain
discount_factor <- 1/((1+r)^c(1:t))
# Make a numeric matrix of the discount_factor with the dimensions of
# rows = 6425 (number of sites in the domain), columns = 10 (time horizon)
r_iy_aqua <- do.call(rbind, replicate(length(fulldomain), discount_factor, simplify=FALSE))

# Calculate the 10-year Net Present Value (NPV) and annuities for each form of aquaculture
# and halibut (the only impacted sector with direct monetary value)
  # Function to calculate the NPV/annuities for Mussel and Kelp
  Value.MK <- function(yield,upfront.cost,annual.cost,price){
    revenue <- do.call(cbind, replicate(t, yield * price, simplify=FALSE))
    cost <- cbind(upfront.cost + annual.cost,
      do.call(cbind, replicate(t - 1, annual.cost, simplify=FALSE)))
    profit <- (revenue - cost) * r_iy_aqua
    profit[profit < 0] <- 0
    NPV <- apply(profit, FUN = sum, MARGIN = 1)
    Annuity = (r*NPV)/(1-((1+r)^-t))
    return(list(NPV = NPV,Annuity = Annuity))
  }
  # Function to calculate the NPV/annuities for Finfish
  Value.F <- function(yield,costs,price){
    revenue <- do.call(cbind, replicate(t, yield * price, simplify=FALSE))
    cost <- sector_data.df$fish.annual.operating.costs
    profit <- (revenue - cost) * r_iy_aqua
    profit[profit < 0] <- 0
    NPV <- apply(profit, FUN = sum, MARGIN = 1)
    Annuity = (r*NPV)/(1-((1+r)^-t))
    return(list(NPV = NPV,Annuity = Annuity))
  }
# Calculate respective NPV/annuities
# Mussel, fixed price of $3.30 per kg
  M <- Value.MK(sector_data.df$mussel.yield,
    sector_data.df$mussel.upfront.cost,
    sector_data.df$mussel.annual.operating.cost,3.3)
# Finfish, fixed price of $8.00 per kg
  F <- Value.F(sector_data.df$fish.yield,
    sector_data.df$fish.upfront.cost,
    unique(sector_data.df$fish.price[sector_data.df$fish.price>0]))
# Kelp, fixed price of $3.00 per kg
  K <- Value.MK(sector_data.df$kelp.yield,
    sector_data.df$kelp.upfront.cost,
    sector_data.df$kelp.annual.operating.cost,3)
# Remove unprofitable sites and generate seperate vectors for
# Those sites in which ventures will be profitable for at least one type
# of aquaculture (var Aqua.Full.Domain),
Aqua.Full.Domain <- data.frame(M$Annuity,F$Annuity,K$Annuity)
Aqua.Full.Domain.Logical <- apply(1 * (Aqua.Full.Domain > 0),FUN=sum,MARGIN=1) > 0
<<<<<<< HEAD
# Sites that will be profitable for mussel (M.V_n_i_p)
M.V_n_i_p <- M$Annuity[Aqua.Full.Domain.Logical]
# Sites that will be profitable for finfish (F.V_n_i_p)
F.V_n_i_p <- F$Annuity[Aqua.Full.Domain.Logical]
# Sites that will be profitable for kelp (K.V_n_i_p)
K.V_n_i_p <- K$Annuity[Aqua.Full.Domain.Logical]
=======
# Sites that will be profitable for mussel (M.V)
M.V <- M$Annuity[Aqua.Full.Domain.Logical]
# Sites that will be profitable for finfish (F.V)
F.V <- F$Annuity[Aqua.Full.Domain.Logical]
# Sites that will be profitable for kelp (K.V)
K.V <- K$Annuity[Aqua.Full.Domain.Logical]
>>>>>>> e9a6acd68e4a0fd08b18829858208289c9350661

# Run the Halibut fishing model and then load the results
if(readline("Run halibut model or load results Y/N? ") == 'Y'){
  print("Launching MATLAB.....");
  system2('/Applications/MATLAB_R2016b.app/bin/matlab',
    args = c('-nodesktop','-noFigureWindows','-nosplash','-r',
    "run\\(\\'~/MSP_Model/Scripts/Halibut/Tuner_free_params_v4.m\\'\\)"))
  # run_matlab_script(paste0(wkdir,'/MSP_Model/Scripts/Halibut/Tuner_free_params_v4.m'))
}
<<<<<<< HEAD
H <- read_excel(paste0(wkdir,'/MSP_Model/Output/Target_FID_and_Yi_fulldomain_NPV_at_MSY_noAqua.xlsx'))

# Load Viewshed Data
F.Viewshed <- as.numeric(gsub(",", "", sector_data.df$res_views_8k)) + as.numeric(gsub(",", "",sector_data.df$park_views_8k))
MK.Viewshed <- as.numeric(gsub(",", "", sector_data.df$res_views_3k)) + as.numeric(gsub(",", "",sector_data.df$park_view_3k))

# Calculate disease impacts or load previously saved results
# Load pathogen connectivity matrix
# Temp file
filename <- paste(tempfile(tmpdir = paste0(wkdir,'/MSP_Model/Input/Data/')),".mat",sep="")
# Write a .mat file with the filtered connectivity matrix
writeMat(filename,
eig = readMat(paste0(wkdir,
  '/MSP_Model/Input/Data/disease_connect_matrix.mat'))$disease.connect.matrix[F$Annuity > 0,F$Annuity > 0])
# Character vector to send to matlab from R
=======
H.V  <- (r*read.csv(paste0(wkdir,'/MSP_Model/Output/Data/Target_FID_and_Yi_fulldomain_NPV_at_MSY_noAqua.csv'),header=FALSE)[,2][Aqua.Full.Domain.Logical])/(1-((1+r)^-t))

# Load Viewshed Data
V_F.V <- (as.numeric(gsub(",", "", sector_data.df$res_views_8k)) + as.numeric(gsub(",", "",sector_data.df$park_views_8k)))[Aqua.Full.Domain.Logical]
V_MK.V  <- (as.numeric(gsub(",", "", sector_data.df$res_views_3k)) + as.numeric(gsub(",", "",sector_data.df$park_view_3k)))[Aqua.Full.Domain.Logical]

# Load Benthic Data, for cells which are not developable for fish aqua set to NA
B.V <- rep(NA,times = length(F.V))
B.V[F.V > 0] <- sector_data.df$TOC.flux[Aqua.Full.Domain.Logical][F.V > 0]

# Run the eigenvector centrality diseaase model in MATLAB and then load the results.
# Write a .mat file with the filtered connectivity matrix
filename <- paste(tempfile(tmpdir = paste0(wkdir,'/MSP_Model/Input/Data')),".mat",sep="")
writeMat(filename,
eig = readMat(paste0(wkdir,
  '/MSP_Model/Input/Data/disease_connect_matrix.mat'))$disease.connect.matrix[F$Annuity > 0,F$Annuity > 0])
# Character vector to send to MATLAB from R. The function eigencentrality is derived from http://strategic.mit.edu/downloads.php?page=matlab_networks
>>>>>>> e9a6acd68e4a0fd08b18829858208289c9350661
code <- c("cd(strcat(pwd,\'/MSP_Model/Scripts/\'));",paste0('load \'',filename,'\';'),'d = abs(eigencentrality(eig));',
'save(\'tmp.mat\',\'d\')')
# Send arguments to matlab
run_matlab_code(code)
# Read the Mat file and remove the temporary one
<<<<<<< HEAD
D <- readMat(paste0(wkdir,'/MSP_Model/Scripts/tmp.mat'))$d
system2('rm',args = paste0(wkdir,'/MSP_Model/Scripts/tmp.mat'))

















# Calculate R_max for viewshed for each developable site
Vi <- apply(cbind(F.Viewshed,MK.Viewshed),MARGIN = 1, FUN = max)[Aqua.Full.Domain.Logical]
# Find the maximum response across all sites
V.R_max <- max(Vi)
# Invert Finfish responses and Mussel/Kelp responses
VF.V_n_i_p <- V.R_max - F.Viewshed[Aqua.Full.Domain.Logical]
VMK.V_n_i_p <- V.R_max - MK.Viewshed[Aqua.Full.Domain.Logical]
## Load Benthic Data
Bi <- df$TOC.flux[F.NPV > 0]
# Find the maximum response across all sites
B.R_max <- max(Bi)
# Invert Finfish responses and Mussel/Kelp responses
B.V_n_i_p <- rep(0, times = length(which(Aqua.Full.Domain.Logical)))
B.V_n_i_p[F.NPV[Aqua.Full.Domain.Logical] > 0] <- B.R_max - Bi
# Disease
## Load disease connectivity matrix, which will be used as the adjacency matrix
# for the disease propagation network
disease_mat <- readMat(paste0(wkdir,'/MSP_Model/Input/Data/disease_connect_matrix.mat'))$disease.connect.matrix
# Remove all sites which cannot be developed for finfish
disease_mat_finfish <- disease_mat[F$Annuity > 0,F$Annuity > 0]
system2('/Applications/MATLAB_R2016b.app/bin/matlab',
  args = c('-nodesktop','-noFigureWindows','-nodisplay','-nosplash',
  '-r \\"try, r'))
graph <- graph_from_adjacency_matrix(disease_mat_finfish, weighted = T, mode = 'undirected')
Di_compare <- eigen_centrality(graph, directed = F, weights = E(graph)$weight)$vector
# plot(1:392,Di_compare/sum(Di_compare))
# points(1:392,Di/sum(Di),col='red')
Di <- read.csv(file = paste0(model_directory,'Data/Raw_Patch_Data.csv'))[F.NPV[Aqua.Full.Domain.Logical] > 0,7]
# Find the maximum response across all sites
D.R_max <- max(Di)
# Invert Finfish responses and Mussel/Kelp responses
D.V_n_i_p <- rep(0, times = length(which(Aqua.Full.Domain.Logical)))
D.V_n_i_p[F.NPV[Aqua.Full.Domain.Logical] > 0] <- D.R_max - Di
## Load Benthic Data
## Disease Model

## Invert impacted sector values
## Scale annuities for all sectors
## Run MSP Model
## Run Dynamic Halibut Model
## Plot Results









# run_matlab_script(paste0(wkdir,'/MSP_Model/Scripts/Halibut_Model_2016/Halibut_tuner_free_params_v4.m'))
#
# # Insert the directory in which the model folder is located
# model_directory <- '~/Desktop/Aquaculture_MSP_Model/'
# # Insert directory in which MATLAB is located
# matlab_directory <- '/Applications/MATLAB_R2016b.app/bin'
#
# # Run MSP tradeoff model in MATLAB
# run_matlab_script('~/Desktop/Aquaculture_Paper/Code/TOA_AquaMSP_CrowCode_v1NaN.m')
# run_matlab_script('~/Desktop/Aquaculture_Paper/Code/EFpayoffs_AquaMSP_CrowCode_v2.m')
#
# n.sector <- 7 # Number of sectors
# epsilon <- 0.2 # Stepsize of sector weights
#
# t <- 10 # Time Horizon
# r <- 0.05 # Discount rate
#
# # Load sector data calculated from sector-specific bioeconomic models
#   df <- read.csv(file=paste0(model_directory,'Data/SeaGrant_data_complete_2015.csv'),stringsAsFactors=FALSE)
#   fulldomain <- df$TARGET_FID # Model domain
#   discount_factor <- 1/((1+r)^c(1:t))
#   r_iy_aqua <- do.call(rbind, replicate(length(fulldomain), discount_factor, simplify=FALSE))
# #### Calculate Responses (R) to Development Options
# ## Sectors which higher responses increase value
#   # Aquaculture
#     Value.KM <- function(yield,upfront.cost,annual.cost,price){
#       revenue <- do.call(cbind, replicate(t, yield * price, simplify=FALSE))
#       cost <- cbind(upfront.cost + annual.cost,
#         do.call(cbind, replicate(t - 1, annual.cost, simplify=FALSE)))
#       profit <- (revenue - cost) * r_iy_aqua
#       profit[profit < 0] <- 0
#       NPV <- apply(profit, FUN = sum, MARGIN = 1)
#       Annuity = (r*NPV)/(1-((1+r)^-t))
#       return(list(NPV = NPV,Annuity = Annuity))
#     }
#     Value.F <- function(yield,costs,price){
#       revenue <- do.call(cbind, replicate(t, yield * price, simplify=FALSE))
#       cost <- df$fish.annual.operating.costs
#       profit <- (revenue - cost) * r_iy_aqua
#       profit[profit < 0] <- 0
#       NPV <- apply(profit, FUN = sum, MARGIN = 1)
#       Annuity = (r*NPV)/(1-((1+r)^-t))
#       return(list(NPV = NPV,Annuity = Annuity))
#     }
#
#     M <- Value.KM(df$mussel.yield,
#       df$mussel.upfront.cost,
#       df$mussel.annual.operating.cost,3.3)
#     F <- Value.F(df$fish.yield,
#       df$fish.upfront.cost,
#       unique(df$fish.price[df$fish.price>0]))
#     K <- Value.KM(df$kelp.yield,
#       df$kelp.upfront.cost,
#       df$kelp.annual.operating.cost,3)
#
#     # Remove unprofitable cells to give the raw
#     Aqua.Full.Domain <- data.frame(M$Annuity,F$Annuity,K$Annuity)
#     Aqua.Full.Domain.Logical <- apply(1 * (Aqua.Full.Domain > 0),FUN=sum,MARGIN=1) > 0
#     # Aqua <- Aqua.Full.Domain[Aqua.Full.Domain.Logical,]
#     M.V_n_i_p <- M$Annuity[Aqua.Full.Domain.Logical]
#     F.V_n_i_p <- F$Annuity[Aqua.Full.Domain.Logical]
#     K.V_n_i_p <- K$Annuity[Aqua.Full.Domain.Logical]
#   # Halibut
#   Halibut_10yrNPVi <- read.csv(file=paste0(model_directory,'Data/Target_FID_and_Yi_fulldomain_NPV_at_MSY_noAqua.csv'));
#   names(Halibut_10yrNPVi) <- c('FID','Yi')
#   H.NPV <- Halibut_10yrNPVi[Aqua.Full.Domain.Logical,2]
#   H.V_n_i_p <- H.NPV
# ## Sectors which higher responses decreases value
#   # Viewshed
#     F.Viewshed <- as.numeric(gsub(",", "", df$res_views_8k)) + as.numeric(gsub(",", "",df$park_views_8k))
#     MK.Viewshed <- as.numeric(gsub(",", "", df$res_views_3k)) + as.numeric(gsub(",", "",df$park_view_3k))
#     # Calculate R_max for viewshed for each developable site
#     Vi <- apply(cbind(F.Viewshed,MK.Viewshed),MARGIN = 1, FUN = max)[Aqua.Full.Domain.Logical]
#     # Find the maximum response across all sites
#     V.R_max <- max(Vi)
#     # Invert Finfish responses and Mussel/Kelp responses
#     VF.V_n_i_p <- V.R_max - F.Viewshed[Aqua.Full.Domain.Logical]
#     VMK.V_n_i_p <- V.R_max - MK.Viewshed[Aqua.Full.Domain.Logical]
#   # Benthic
#     Bi <- df$TOC.flux[F.NPV > 0]
#     # Find the maximum response across all sites
#     B.R_max <- max(Bi)
#     # Invert Finfish responses and Mussel/Kelp responses
#     B.V_n_i_p <- rep(0, times = length(which(Aqua.Full.Domain.Logical)))
#     B.V_n_i_p[F.NPV[Aqua.Full.Domain.Logical] > 0] <- B.R_max - Bi
#   # Disease
#     ## Load disease connectivity matrix, which will be used as the adjacency matrix
#     # for the disease propagation network
#     disease_mat <- readMat(paste0(model_directory,'Data/disease_connect_matrix.mat'))$disease.connect.matrix
#     # Remove all sites which cannot be developed for finfish
#
#     ## 12/29/16 differences between Matlab and R eigenvector centrality metrics
#     disease_mat_finfish <- disease_mat[F$Annuity > 0,F$Annuity > 0]
#     graph <- graph_from_adjacency_matrix(disease_mat_finfish,diag = T, weighted = T, mode = 'undirected')
#     Di_compare <- eigen_centrality(graph, directed = T, weights = E(graph)$weight)$vector
#     # plot(1:392,Di_compare/sum(Di_compare))
#     # points(1:392,Di/sum(Di),col='red')
#     Di <- read.csv(file = paste0(model_directory,'Data/Raw_Patch_Data.csv'))[F.NPV[Aqua.Full.Domain.Logical] > 0,7]
#     # Find the maximum response across all sites
#     D.R_max <- max(Di)
#     # Invert Finfish responses and Mussel/Kelp responses
#     D.V_n_i_p <- rep(0, times = length(which(Aqua.Full.Domain.Logical)))
#     D.V_n_i_p[F.NPV[Aqua.Full.Domain.Logical] > 0] <- D.R_max - Di
# ### Scale Sectors
#   Scaling <- function(Sector){Scaled_Sector <- Sector / max(Sector)}
#    M.X <- Scaling(M.V_n_i_p)
#    F.X <- Scaling(F.V_n_i_p)
#    K.X <- Scaling(K.V_n_i_p)
#    H.X <- Scaling(H.V_n_i_p)
#    VF.X <- Scaling(VF.V_n_i_p)
#    VMK.X <- Scaling(VMK.V_n_i_p)
#    D.X <- Scaling(D.V_n_i_p)
#
#   Scaled_Data <- data.frame(M_X = M.X, F_X = F.X, K_X = K.X,
#     H_X = H.X, VF_X = VF.X, VMK_X = VMK.X, D_X = D.X)
# ## Analyze data to make extract all values
# # Load old version results
# raw_values.df <- read.csv('~/Desktop/CrowTOv1/Raw_Patch_Data.csv')
# names(raw_values.df) <- c('M','F','K','H','V_F','V_MK','B','D')
# write.csv(file = '~/Desktop/Aquaculture_Paper/Code/Raw_Patch_Data.csv',raw.annuities)
#
# data.df <- read.csv(file='~/Desktop/Code/MSP Planning Results April 2016/Dynamic_Values_Export.csv',header=F)
# names(data.df) <- c('M','F','K','H','V','B','D')
# plans.df <- read.csv(file='~/Desktop/Code/MSP Planning Results April 2016/Static_plans.csv',header=F)
#
# # Conventional Models
# Unconstrained_data <- read.csv(file='~/Desktop/Code/MSP Planning Results April 2016/Unconstrained_Dynamic_Values_April.csv',header=F)
# names(Unconstrained_data)=c('Mussel','Finfish','Kelp','Halibut','Viewshed','Benthic','Disease')
# Constrained_data <- read.csv(file='~/Desktop/Code/MSP Planning Results April 2016/Constrained_Dynamic_Values_April.csv',header=F)
# names(Constrained_data)=c('Mussel','Finfish','Kelp','Halibut','Viewshed','Benthic','Disease')
#
# monetry_values.df <- data.frame(M = data.df$M * sector_totals.df[1],
#   F = data.df$F * sector_totals.df[2],
#   K = data.df$K * sector_totals.df[3],
#   H = data.df$H * sector_totals.df[4])
# # Conversion of NPV to equivalent annuity
# # Calculate annuities
# annuities <- monetry_values.df %>% mutate(M = (r*M)/(1-((1+r)^-t)),
#                 F = (r*F)/(1-((1+r)^-t)),
#                   K = (r*K)/(1-((1+r)^-t)),
#                   H = (r*H)/(1-((1+r)^-t)))
# # Calculate the annuitiy for > 25% of Mussel and loss of annuity of less than 1%
# data.df %>% filter(M > .25, 1 - H < .01) %>% select(M,H) %>%
#       mutate(M = M * sector_totals.df[1],H = (1-H) * sector_totals.df[4])
#
# # Calculate the case study
# Case_Study.df <- data.df %>% mutate(ID = 1:nrow(data.df)) %>% filter(M >= .05, F >= .05, K >= .05, H >= .95, V >= .95, B >= .95, D >= .95)
# Case_Study.df[576,] %>% select(M,F,K) %>% mutate(M = M * sector_totals.df[1], F = F * sector_totals.df[2], K = K * sector_totals.df[3])
# summary(factor(plans.df[,Case_Study.df$ID[576]]))
#
# # Calculate sector differences
# color.vector=c(rep('coral',length.out=nrow(data.df)),rep('purple',length.out=1061),rep('green',length.out=1061))
#  panel.EF<-function(x, y, itor=0, epsilon=.001, bg = NA, pch = 20, cex = .01, ...)
#   {
#     x.MSP=x[color.vector=='coral']
#     y.MSP=y[color.vector=='coral']
#
#     # points(x.MSP,y.MSP, pch = 16, col = alpha("lightblue1",1/100),cex = cex)
#     x.U=x[color.vector=='purple']
#     y.U=y[color.vector=='purple']
#     x.S=x[color.vector=='green']
#     y.S=y[color.vector=='green']
#     x.EF=NULL
#     y.EF=NULL
#     alpha.mat.tmp=seq(from=0,by=epsilon,to=1)
#     # MSP
#     for(itor in 1:length(alpha.mat.tmp)){
#       alpha.tmp=alpha.mat.tmp[itor]
#       A=(alpha.tmp*x.MSP)+((1-alpha.tmp)*y.MSP)
#       I=which(A==max(A))
#       x.EF[itor]=max(unique(x.MSP[I]))
#       I.tmp.x=which(x.MSP==max(unique(x.MSP[I])))
#       I.tmp.y=which(y.MSP[I.tmp.x]==max(unique(y.MSP[I.tmp.x])))
#       y.EF[itor]=unique(y.MSP[I.tmp.x[I.tmp.y]])}
#     x.EF.original=x.EF;y.EF.original=y.EF;
#     if(length(unique(x.EF.original))!=1&length(unique(x.EF.original))!=1){
#       EF.inter=approx(x.EF.original,y=y.EF.original,n=length(alpha.mat.tmp))
#       x.EF=EF.inter$x;y.EF=EF.inter$y;
#     }else{
#     }
#     # lines(sort(x.EF),y.EF[order(x.EF)],col="midnightblue",lwd=2,lty=1)
#     # lines(sort(x.U),y.U[order(x.U)],col = "mediumorchid1",lwd=2,lty=1)
#     # lines(sort(x.S),y.S[order(x.S)],col = "coral1",lwd=2,lty=1)
#     return(cbind(x.EF,y.EF))
#   }
#
# M_H.EF <- data.frame(panel.EF(data.df$M,data.df$H)) %>%
#     mutate(M_raw.EF = x.EF * sector_totals.df[1],
#     H_raw.EF = y.EF * sector_totals.df[4])
# M_H.S <- data.frame(Constrained_data[,1], Constrained_data[,4], Constrained_data[,1] * sector_totals.df[1], Constrained_data[,4] * sector_totals.df[4])
# names(M_H.S) <- c('M.percentage', 'H.percentage', 'M.annuity', 'H.annuity')
#
# M_H.EF %>% filter(x.EF == .440)
# tail(M_H.S)
#
#
# # Load Scenario Data
# # print('Loading Data........')
# # V1_JS.data <- read.csv(file = '~/Desktop/Code/V1_Policy_JS.csv',header = F)
# # V2_JS.data <- read.csv(file = '~/Desktop/Code/V2_Policy_JS.csv',header = F)
# # V2_CW.data <- read.csv(file = '~/Desktop/Code/V2_Policy_CW.csv',header = F)
#
# # freq_scenario <- function(data){
# #   List.Data <- list()
# #   fx.t <- proc.t()
# #   print('This may take a while...........')
# #   for(index in 1:ncol(data)){
# #     if(index %in% seq(from = 0, to = ncol(data), by = 10000)){
# #       print(index)
# #       ptm <- proc.t()
# #     }
# #     tab.tmp <- data.frame(table(factor(data[,index],
# #       levels = c(1:4),
# #       labels = c('No Development','Mussel','Finfish','Kelp'))))
# #     names(tab.tmp) <- c('Policy','Frequency')
# #     list.name <- paste("Scenario_", index, sep="")
# #     List.Data[[list.name]] <- tab.tmp %>% mutate(Scenario = rep(index,ts = length(tab.tmp)))
# #   }
# #   new.data <- rbindlist(List.Data)
# #   if(index %in% seq(from = 0, to = ncol(data), by = 10000)){
# #      return(print(proc.t()-ptm))
# #   }
# #   print(paste('t elapsed',proc.t() - fx.t))
# #   return(new.data)
# # }
# # print('Arranging Data........')
# # SW_Version.df <- freq_scenario(V1_JS.data) %>% mutate(Version = 'Old Approach')
# # JS_Version.df <- freq_scenario(V2_JS.data) %>% mutate(Version = 'New Approach JS')
# # CW_Version.df <- freq_scenario(V2_CW.data) %>% mutate(Version = 'New Approach CW')
#
# # # Combine all
# # Version.df <- bind_rows(SW_Version.df %>% mutate(Version = 'Old Approach'),
# #     JS_Version.df %>% mutate(Version = 'New Approach JS'),
# #     CW_Version.df %>% mutate(Version = 'New Approach CW'))
#
# # head(Version.df %>% filter(Version == 'Old Approach'))
# # head(Version.df %>% filter(Version == 'New Approach JS'))
# # head(Version.df %>% filter(Version == 'New Approach CW'))
# # # Plots to compare percentage
# # png(filename = '~/Desktop/Code/Compare.png')
# # Version.df %>% ggplot(aes(x = Scenario, y = Frequency)) +
# #   geom_line(aes(colour = Policy)) +
# #   facet_grid(Version~.) +
# #   theme_minimal()
# # dev.off()
#
# # write.csv(file = '~/Desktop/Code/Version.csv',Version.df)
#
#
#
#
#
# # ## Figures
# # This is the most concise figure script, created 12/28/16 by JS
# # Makes all five primary figures for MSP Aquaculture Paper
# # Uses maps created by Becca Gentry in GIS
# current.directory.data <- paste0(model_directory,'Data')
# ## Figure dimensions
# width=7
# height=5.5
# res=2400
# units='in'
# ## Set Directory
# cols <- c("Mussel"="Blue","Finfish"="Salmon","Kelp"="Green","Halibut"="Burlywood","Viewshed"='Cyan',"Benthic"='Orange',"Disease"='Black')
# text.size<-12
# patch.size=1.5
# aMatrix <- read.csv(file='~/Desktop/Code/MSP Planning Results April 2016/aMatrix.csv',header=F)
# # Static.plans.data <- read.csv(file='Static_plans.csv',header=F)
# Static.values.data <- read.csv(file='~/Desktop/Code/MSP Planning Results April 2016/Dynamic_Values_Export.csv',header=F)
# colnames(Static.values.data)=c('Mussel','Finfish','Kelp','Halibut','Viewshed','Benthic','Disease')
# # Static.plans.case.study.data <- read.csv(file='~/Desktop/Code/MSP Planning Results April 2016/April Results/Static_plans_case_study.csv',header=F)
# # Static.values.case.study.data <- read.csv(file='~/Desktop/Code/MSP Planning Results April 2016/Static_values_case_study.csv',header=F)
# # Static.percentage.case.study.data <- read.csv('~/Desktop/Code/MSP Planning Results April 2016/Static_percentage_case_study.csv',header=F)
# # colnames(Static.values.data)=c('Mussel','Finfish','Kelp','Halibut','Viewshed','Benthic','Disease')
# # colnames(Static.values.case.study.data)=c('Mussel','Finfish','Kelp','Halibut','Viewshed','Benthic','Disease')
# # data.MSP = list(aMatrix = aMatrix,Static.plans = Static.plans.data,Static.plans.case.study = Static.plans.case.study.data,
# #               Static.values.case.study = Static.values.case.study.data,Static.percentage.case.study = Static.percentage.case.study.data)
# # Case Study Stuff
# I<-which(Static.values.data$Mussel>=.05&Static.values.data$Finfish>=.05&Static.values.data$Kelp>=.05&
#                            Static.values.data$Halibut>=.95&Static.values.data$Viewshed>=.95&Static.values.data$Benthic>=.95&
#                            Static.values.data$Disease>=.95)
# df <- read.csv(file='~/Desktop/Code/SeaGrant_data_complete_2015.csv',stringsAsFactors=FALSE)
# V1 <- read.csv(file='~/Desktop/Code/V1.csv',header=FALSE)
# names(df)
# paste('Cheapest Mussel Farm Costs',df %>% filter(V1 == 1) %>%
#   select(mussel.annual.operating.costs) %>%
#   filter(mussel.annual.operating.costs > 0) %>% arrange(mussel.annual.operating.costs) %>% filter(mussel.annual.operating.costs == min(mussel.annual.operating.costs)),'annually')
#
# # Case.Study.August <- read.csv(file='Case Study August.csv',header=F)
# # Static.values.data[I[576],]
# # summary(factor(Static.plans.data[,I[576]]))
# # Conventional Models
# Unconstrained_data <- read.csv(file='~/Desktop/Code/MSP Planning Results April 2016/Unconstrained_Dynamic_Values_April.csv',header=F)
# names(Unconstrained_data)=c('Mussel','Finfish','Kelp','Halibut','Viewshed','Benthic','Disease')
# Constrained_data <- read.csv(file='~/Desktop/Code/MSP Planning Results April 2016/Constrained_Dynamic_Values_April.csv',header=F)
# names(Constrained_data)=c('Mussel','Finfish','Kelp','Halibut','Viewshed','Benthic','Disease')
#
# # Unconstrained_data <- read.csv(file='Unconstrained_Static_Values_April.csv',header=F)
# # names(Unconstrained_data)=c('Mussel','Finfish','Kelp','Halibut','Viewshed','Benthic','Disease')
# # Constrained_data <- read.csv(file='Constrained_Static_Values_April.csv',header=F)
# # names(Constrained_data)=c('Mussel','Finfish','Kelp','Halibut','Viewshed','Benthic','Disease')
#
# # Pure Profit Conventional Models
# # Unconstrained.Pure.Profit.data <- read.csv(file='Pure_Profit_Unconstrained_Dynamic_Values_April.csv',header=F)
# # names(Unconstrained.Pure.Profit.data)=c('Mussel','Finfish','Kelp','Halibut','Viewshed','Benthic','Disease')
# # Constrained.Pure.Profit.data <- read.csv(file='Pure_Profit_Constrained_Dynamic_Values_April.csv',header=F)
# # names(Constrained.Pure.Profit.data)=c('Mussel','Finfish','Kelp','Halibut','Viewshed','Benthic','Disease')
# Master.matrix.max=rbind(Static.values.data,Unconstrained_data,Constrained_data)
# names(Master.matrix.max)<-c('Mussel','Finfish','Kelp','Halibut','Viewshed','Benthic','Disease')
# color.vector.max=c(rep('coral',length.out=nrow(Static.values.data)),rep('purple',ts=nrow(Unconstrained_data)),rep('green',ts=nrow(Constrained_data)))
#
# # Calculate new annuity values
# raw_values.df <- read.csv('~/Desktop/CrowTOv1/Raw_Patch_Data.csv')
# names(raw_values.df) <- c('M','F','K','H','V_F','V_MK','B','D')
# r <- 0.05
# t <- 10
# sector_totals.df <- apply(raw_values.df %>% select(M,F,K,H) %>%
#   mutate(M = (r*M)/(1-((1+r)^-t)),
#                 F = (r*F)/(1-((1+r)^-t)),
#                   K = (r*K)/(1-((1+r)^-t)),
#                   H = (r*H)/(1-((1+r)^-t))),MARGIN=2,FUN = sum)
# raw.annuities <- raw_values.df %>% select(M,F,K,H,V_F,V_MK,B,D) %>%
#   mutate(M = (r*M)/(1-((1+r)^-t)),
#                 F = (r*F)/(1-((1+r)^-t)),
#                   K = (r*K)/(1-((1+r)^-t)),
#                   H = (r*H)/(1-((1+r)^-t)))
# M_group <- gsub(paste(c('\\(','\\]'),collapse = '|'),'',
#   names(split(raw.annuities$M,cut(raw.annuities$M,seq(min(raw.annuities$M[which(raw.annuities$M > 0)]),
#     max(raw.annuities$M),length.out=9)))))
# F_group <- gsub(paste(c('\\(','\\]'),collapse = '|'),'',
#   names(split(raw.annuities$F,cut(raw.annuities$F,seq(min(raw.annuities$F[which(raw.annuities$F > 0)]),
#     max(raw.annuities$F),length.out=9)))))
# K_group <- gsub(paste(c('\\(','\\]'),collapse = '|'),'',
#   names(split(raw.annuities$K,cut(raw.annuities$K,seq(min(raw.annuities$K[which(raw.annuities$K > 0)]),
#     max(raw.annuities$K),length.out=9)))))
# #### Plot Data
# current.directory.scripts='~/Desktop/Code/MS Figures/'
# setwd(current.directory.scripts)
# theme = theme(plot.margin = unit(c(.2,.2,.2,.2), units = "lines"),
#                 axis.text = element_blank(),
#                 axis.title = element_blank(),
#                 axis.ticks = element_blank(),
#                 axis.ticks.length = unit(0, "lines"),
#                 axis.ticks.margin = unit(0, "lines"),
#                 panel.background=element_rect(fill="white"),
#                 panel.grid=element_blank(),
#                 plot.title=element_text(hjust=0))
#     labs = labs(x = NULL, y = NULL)
# # png(paste0(current.directory.scripts,'MS_Figures.png'),onefile = T,units=units,width=width, height=height, res=res)
# for(itor in 1:2){
# if(itor==1){
#   png(paste0(current.directory.scripts,'PNGs/Fig 1.png'),units=units,width=width, height=height, res=res)
# }else{
#   pdf(paste0(current.directory.scripts,'PDFs/Fig 1.pdf'), width=width, height=height,paper='legal')
# }
#   # img <- readTIFF("fig1_Stevens_v3.tif",native=T,info=T)
#   # g <- rasterGrob(img, interpolate=TRUE)
#   img <- readPNG("~/Desktop/Code/Fig1A Capture.png",native=T,info=T)
#   g <- rasterGrob(img, interpolate=TRUE)
#
#   a<-qplot(1:10, 1:10, geom="blank") + annotation_custom(g, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) + ggtitle('A') +
#     theme + labs
#
#   img_mussel <- readPNG("~/Desktop/Code/MusselValueApril.png",native=T,info=T)
#   g_mussel <- rasterGrob(img_mussel, interpolate=TRUE)
#
#   foo <- 5.58
#   store <- NULL
#   for(itor in 1:9){
#     foo <- foo - .341
#     store[itor] <- foo
#     # print(a)
#   }
#
#   b <- qplot(1:10, 1:10, geom="blank") + ggtitle('B') + annotation_custom(g_mussel, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +
#     theme +
#     # annotate(geom='text',x = c(1.75,1.83,1.83,1.83,1.92,1.92,1.92,1.92,1.89),
#     #   y = store,label=c('0-6','6-6.5','6.5-7','7-7.5','7.5-8.0','8.0-8.5','8.5-9.0','9.0-9.5','9.5-10')) +
#     labs
#   # ggsave(filename = "~/Desktop/Code/MusselValueApril_Edit.png",b,device = 'png',dpi = res)
#
#   foo <- 5.57
#   store <- NULL
#   for(itor in 1:9){
#     foo <- foo - .35
#     store[itor] <- foo
#     # print(a)
#   }
#
#   img_finfish <- readPNG("~/Desktop/Code/FishValueApril.png",native=T,info=T)
#   g_finfish <- rasterGrob(img_finfish, interpolate=TRUE,just='center')
#
#   c<-qplot(1:10, 1:10, geom="blank") + ggtitle('C') + annotation_custom(g_finfish, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +
#     theme +
#     # annotate(geom='text',x = c(1.87,1.97,1.92,1.92,1.92,1.92,1.92,1.92,1.92),
#     #   y = store,label=c('0-0.01','0.01-0.4','0.4-0.8','0.8-1.2','1.2-1.6','1.6-2.0','2.0-2.4','2.4-2.8','2.8-3.2')) +
#     labs
#   # ggsave(filename = "~/Desktop/Code/FishValueApril_Edit.png",c,device = 'png',dpi = res)
#
#   img_kelp <- readPNG("~/Desktop/Code/KelpValueApril.png",native=T,info=T)
#   g_kelp <- rasterGrob(img_kelp, interpolate=TRUE, just='center')
#
#   # foo <- 5
#   # store <- NULL
#   # for(itor in 1:9){
#   #   foo <- foo - 0.375
#   #   store[itor] <- foo
#   #   # print(a)
#   # }
#
#   foo <- 5.58
#   store <- NULL
#   for(itor in 1:10){
#     foo <- foo - .315
#     store[itor] <- foo
#     # print(a)
#   }
# # annotate(geom='text',x = c(rep(.75,ts=7),.75+.09,.75+.09)
# # store[length(store) - 2] <- store[length(store) - 2] + .1
# # store[length(store) - 2] <- store[length(store) - 1] + .2
#   d<-qplot(1:10, 1:10, geom="blank") + ggtitle('D') + annotation_custom(g_kelp, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +
#     theme +
#     # annotate(geom='text',x = c(1.75,1.83,1.83,1.83,1.92,1.92,1.92,1.92,1.89),
#     #   y = store,label=c('0-3','3-4','4-5','5-6','6-7','7-8','8-9','9-10','10-12'),size = 3.85) +
#      # annotate(geom='text',x = c(1.75,1.75,1.75,1.75,1.75,1.75,1.75,1.75,1.75,1.80),
#      #  y = store,label=c('0-1','1-2','2-3','3-4','4-5','5-6','6-7','7-8','8-9','10-12')) +
#     labs #+
#     # scale_x_continuous(breaks = seq(0, 10, .25)) +
#     # scale_y_continuous(breaks = seq(0, 10, .25)) +
#     # theme(panel.ontop=TRUE,panel.background = element_rect(colour = NA,fill="transparent"),panel.grid.minor = element_blank(),panel.grid.major=element_line(color = 'black'))
#   # ggsave(filename = "~/Desktop/Code/KelpValueApril_Edit.png",d)
#
#   grid.arrange(a,b,c,d,ncol=2,nrow=2)
#                  # bottom = textGrob(expression(bold('Fig 1: ')~plain('Study domain, spatial constraints and potential value for aquaculture development. (A) Spatial constraints to aquaculture development in the Southern California Bight. (B-D) Potential value (10- year NPV) in each developable cell for mussel, finfish, and kelp aquaculture sectors.')),
#                  #                   x=1,just='left'))
#   if(itor == 1){
#       dev.off()
#     }else{
#       dev.off()
#     }
# }
# # Figure 2
# # pdf(paste0(current.directory.scripts,'Fig 2.pdf'),width=8, height=6.4,paper='legal')
# # img_MSP <- readPNG(paste0(current.directory.scripts,'Fig 2.png'),native=T,info=T)
# #   g_MSP <- rasterGrob(img_MSP, interpolate=TRUE, just='center')
# #   qplot(1:10, 1:10, geom="blank") + annotation_custom(g_MSP, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +
# #     theme + labs
# # dev.off()
#   # png(paste0(current.directory.scripts,'MS_Figures.png'),onefile = T, width=8, height=6.4,res=res,units=units)
# for(itor in 1){
#   if(itor==1){
#     png(paste0(current.directory.scripts,'PNGs/Fig 2.png'),width=8, height=6.4,res=res,units=units)
#   }else{
#     pdf(paste0(current.directory.scripts,'PDFs/Fig 2.pdf'),width=8, height=6.4,paper='legal')
#   }
#   panel.EF<-function (x, y, itor=0, epsilon=.001, bg = NA, pch = 20, cex = .01, ...)
#   {
#     x.MSP=x[color.vector=='coral']
#     y.MSP=y[color.vector=='coral']
#
#     points(x.MSP,y.MSP, pch = 16, col = alpha("dodgerblue",1/75),cex = cex/2)
#     x.U=x[color.vector=='purple']
#     y.U=y[color.vector=='purple']
#     x.S=x[color.vector=='green']
#     y.S=y[color.vector=='green']
#     x.EF=NULL
#     y.EF=NULL
#     alpha.mat.tmp=seq(from=0,by=epsilon,to=1)
#     # MSP
#     for(itor in 1:length(alpha.mat.tmp)){
#       alpha.tmp=alpha.mat.tmp[itor]
#       A=(alpha.tmp*x.MSP)+((1-alpha.tmp)*y.MSP)
#       I=which(A==max(A))
#       x.EF[itor]=max(unique(x.MSP[I]))
#       I.tmp.x=which(x.MSP==max(unique(x.MSP[I])))
#       I.tmp.y=which(y.MSP[I.tmp.x]==max(unique(y.MSP[I.tmp.x])))
#       y.EF[itor]=unique(y.MSP[I.tmp.x[I.tmp.y]])}
#     x.EF.original=x.EF;y.EF.original=y.EF;
#     if(length(unique(x.EF.original))!=1&length(unique(x.EF.original))!=1){
#       EF.inter=approx(x.EF.original,y=y.EF.original,n=length(alpha.mat.tmp))
#       x.EF=EF.inter$x;y.EF=EF.inter$y;
#     }else{
#     }
#     lines(sort(x.EF),y.EF[order(x.EF)],col="midnightblue",lwd=2,lty=1)
#     lines(sort(x.U),y.U[order(x.U)],col = "mediumorchid1",lwd=2,lty=1)
#     lines(sort(x.S),y.S[order(x.S)],col = "coral1",lwd=2,lty=1)
#   }
#   #   pdf.options(width = 8, height = 6.4)
#   source('~/Desktop/Code/pairs2.R')
#   color.vector=color.vector.max
#    # Color Vector For Seperating the MSP from Conventional Solutions
#   # sample <- rbind(MM_test.df %>% filter(Set == 'MSP') %>% sample_n(size = 1000),
#   #   MM_test.df %>% filter(Set == 'U') %>% sample_n(size = 500),
#   #   MM_test.df %>% filter(Set == 'C') %>% sample_n(size = 500))
#   pairs2(100*Master.matrix.max,lower.panel=panel.EF,
#          upper.panel=NULL,col=color.vector,cex=0.8,xlim=c(0,100),
#          ylim=c(0,100),pch=16,font.labels=3,cex.axis=1,las=1,xaxp=c(0,100,4),yaxp=c(0,100,2),
#          gap=1)
#   # title(xlab='% of Maximum',line = 1)
#   title(ylab='% of Maximum')
#   par(xpd=T)
#   l1<-legend(.33,1,
#              legend=c('7D Frontier','2D Frontier'),fill=c("lightblue1","midnightblue"),
#              cex=.75,title=expression(bold('Marine Spatial Planning (MSP)')),
#              title.adj = 0, bty = 'n', adj = 0, text.width=.25)
#   l2<-legend(x = l1$rect$left+.0020, y = with(l1$rect, top - h)-.005,
#              legend=c('Constrained','Unconstrained'),fill=c("coral1","mediumorchid1"),
#              cex=.75,title=expression(bold('Conventional Planning ')),
#              title.adj = 0, bty = 'n', adj = 0, text.width=.25)
#   inset.figure.proportion = 1/3
#   inset.figure.dims = c(rep(width*(inset.figure.proportion),ts = 2))
#   subplot(source(file='~/Desktop/Code/Tradeoff Cartoon.R'),x='topright',size = inset.figure.dims, type='plt', par = list(cex.main=2.5, cex = .45, lwd = 1))
#   par(oma=c(0,2,2,0))
#   title('A', adj = 0, outer = T, cex = .75)
#   title(xlab='% of Maximum',line = 3.5)
#   if(itor == 1){
#     dev.off()
#   }else{
#     dev.off()
#   }
# }
# # Figure 3
# for(itor in 1:2){
#   if(itor==1){
#     png(paste0(current.directory.scripts,'PNGs/Fig 3.png'),width=width, height=height,res=res,units=units)
#   }else{
#     pdf(paste0(current.directory.scripts,'PDFs/Fig 3.pdf'),width=width, height=height,paper='legal')
#   }
# # png(paste0(current.directory.scripts,'MS_Figures.png'),onefile = T, width=width, height=height,res=res,units=units)
#
#   img_hot_all <- readPNG("~/Desktop/Code/HotSpots_ALL_Joel.png",native=T,info=T)
#   g_hot_all <- rasterGrob(img_hot_all, interpolate=TRUE)
#
#   a<-qplot(1:10, 1:10, geom="blank") + ggtitle('A') +
#     annotation_custom(g_hot_all, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)+
#     theme + labs
#
#   img_hot_mussel <- readPNG("~/Desktop/Code/HotSpots_Mussels_Joel.png",native=T,info=T)
#   g_hot_mussel <- rasterGrob(img_hot_mussel, interpolate=TRUE)
#
#   b<-qplot(1:10, 1:10, geom="blank") + ggtitle('B') +
#     annotation_custom(g_hot_mussel, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)+
#     theme + labs
#
#   img_hot_finfish <- readPNG("~/Desktop/Code/HotSpots_Fish_Joel.png",native=T,info=T)
#   g_hot_finfish <- rasterGrob(img_hot_finfish, interpolate=TRUE,just='center')
#
#   c<-qplot(1:10, 1:10, geom="blank") + ggtitle('C') +
#     annotation_custom(g_hot_finfish , xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)+
#     theme + labs
#
#   img_hot_kelp <- readPNG("~/Desktop/Code/HotSpots_Kelp_Joel.png",native=T,info=T)
#   g_hot_kelp <- rasterGrob(img_hot_kelp, interpolate=TRUE, just='center')
#
#   d<-qplot(1:10, 1:10, geom="blank") + ggtitle('D') +
#     annotation_custom(g_hot_kelp, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)+
#     theme + labs
#   grid.arrange(a,b,c,d,ncol=2,nrow=2,padding=unit(.1,'line'))
#   dev.off()
# }
# #  pdf(paste0(current.directory.scripts,'Fig 4.pdf'),width=width, height=height,paper='legal')
# # img_MSP <- readPNG(paste0(current.directory.scripts,'PNGs/Fig 4.png'),native=T,info=T)
# # g_MSP <- rasterGrob(img_MSP, interpolate=TRUE, just='center')
# # qplot(1:10, 1:10, geom="blank") + annotation_custom(g_MSP, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +
# #   theme + labs
# # dev.off()
#
# # Figure 4
# for(itor in 1:2){
#   if(itor==1){
#     png(paste0(current.directory.scripts,'PNGs/Fig 4.png'),width=width, height=height,res=res,units=units)
#   }else{
#     pdf(paste0(current.directory.scripts,'PDFs/Fig 4.pdf'),width=width, height=height,paper='legal')
#   }
#   Low.impact.solutions=Static.values.data[I,]
#   LI.names=names(Low.impact.solutions)
#   ID=seq(to=nrow(Low.impact.solutions),from=1,by=1)
#   size.tmp=dim(Low.impact.solutions)
#   Mussel.LI=Low.impact.solutions[,1]
#   Finfish.LI=Low.impact.solutions[,2]
#   Kelp.LI=Low.impact.solutions[,3]
#   Halibut.LI=Low.impact.solutions[,4]
#   View.LI=Low.impact.solutions[,5]
#   Benthic.LI=Low.impact.solutions[,6]
#   Disease.LI=Low.impact.solutions[,7]
#   Sector.LI=c(rep(LI.names[1],ts=size.tmp[1]),
#               rep(LI.names[2],ts=size.tmp[1]),
#               rep(LI.names[3],ts=size.tmp[1]),
#               rep(LI.names[4],ts=size.tmp[1]),
#               rep(LI.names[5],ts=size.tmp[1]),
#               rep(LI.names[6],ts=size.tmp[1]),
#               rep(LI.names[7],ts=size.tmp[1]))
#   for(itor in 1:ncol(Low.impact.solutions)){
#     ID.LI.tmp=1:size.tmp[1]
#     if(itor>1){
#       ID.LI=c(ID.LI,ID.LI.tmp)
#     }else{
#       ID.LI=ID.LI.tmp
#     }
#   }
#   value.LI.tmp=c(Mussel.LI,Finfish.LI,Kelp.LI,Halibut.LI,View.LI,Benthic.LI,Disease.LI)
#   Low.impact.solutions=data.frame(ID.LI,value.LI.tmp,Sector.LI)
#   names(Low.impact.solutions)=c('ID','Value','Sector')
#   Low.impact.solutions$Sector=factor(Low.impact.solutions$Sector, levels=c('Mussel','Finfish','Kelp','Halibut','Viewshed','Benthic','Disease'))
#   p.bar<-ggplot(data = Low.impact.solutions,aes(x=ID,y=Value,fill=Sector,color=Sector))+
#     geom_bar(stat="identity")+facet_grid(.~Sector)+ggtitle('A')+
#     scale_y_continuous(labels = percent,limits=c(0,1))+
#     scale_fill_manual(values=cols)+
#     scale_color_manual(values=cols)+
#     theme(axis.title.x=element_blank(),
#           axis.text.x=element_blank(),
#           axis.ticks.x=element_blank(),
#           panel.background=element_rect(color='white',fill='white'),panel.spacing = unit(.75, "lines"),
#           strip.background=element_rect(fill='white'),strip.text=element_text(size=text.size*.75),
#           panel.grid=element_blank(),axis.title.y=element_text(size=text.size,color="black"),
#           axis.text.y=element_text(size=text.size,color="black"),legend.title=element_text(size=text.size*.75),
#           legend.text=element_text(size=text.size),plot.title=element_text(hjust=0))
#
#   img_case_study_percent <- readPNG("~/Desktop/Code/CaseStudyPercentage_August.png",native=T,info=T)
#   g_case_study_percent <- rasterGrob(img_case_study_percent, interpolate=TRUE,just='center')
#
#   p.map<-qplot(1:10, 1:10, geom="blank") + ggtitle('B') +
#     annotation_custom(g_case_study_percent, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)+
#     theme(plot.margin = unit(c(.2,.2,.4,.2), units = "lines"),axis.text.x = element_blank(),axis.text.y = element_blank(),axis.ticks = element_blank(),
#           axis.title=element_blank(),panel.background=element_rect(fill="white"),panel.grid=element_blank(),
#           plot.title=element_text(hjust=0))
#
#   img_case_study_plan <- readPNG("~/Desktop/Code/CaseStudySpecies_August.png",info=T)
#   g_case_study_plan <- rasterGrob(img_case_study_plan, interpolate=TRUE,just='center')
#
#   p.map.case.study<-qplot(1:10, 1:10, geom="blank") + ggtitle('C') +
#     annotation_custom(g_case_study_plan , xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)+
#     theme(plot.margin = unit(c(.2,.2,.4,.2), units = "lines"),axis.text.x = element_blank(),axis.text.y = element_blank(),axis.ticks = element_blank(),
#           axis.title=element_blank(),panel.background=element_rect(fill="white"),panel.grid=element_blank(),
#           plot.title=element_text(hjust=0))
#
#   p.bar.case.study <- ggplot(data = Low.impact.solutions[Low.impact.solutions$ID==576,],
#                              aes(x=ID,y=Value,fill=Sector,color=Sector))+ggtitle('D')+
#     geom_bar(stat="identity")+facet_grid(.~Sector)+
#     scale_y_continuous(labels = percent,limits=c(0,1))+
#     scale_fill_manual(values=cols)+
#     scale_color_manual(values=cols)+
#     theme(axis.ticks.x=element_blank(),axis.title.x=element_blank(),
#           axis.title.y=element_blank(),axis.text.x=element_blank(),
#           panel.background=element_rect(color='white',fill='white'),
#           strip.background=element_rect(fill='white'),strip.text=element_blank(),
#           panel.grid=element_blank(),
#           axis.text.y=element_blank(),axis.ticks=element_blank(),legend.position='none',plot.title=element_text(hjust=0),plot.margin = unit(c(.2,.2,.2,.2), "cm"))
#   g=ggplotGrob(p.bar.case.study)
#   p.map.case.study.full<-p.map.case.study+annotation_custom(grob=g,xmin=5.5,xmax=10.25,ymin=7.5)
#   grid.arrange(p.bar,p.map,p.map.case.study.full,ncol=2,nrow=2,layout_matrix = rbind(c(1,1),c(2,3)))
#   dev.off()
# }
#   # pdf.options(width = 9.5, height = 7)
# # Figure 5
# for(itor in 1:2){
#   if(itor==1){
#     png(paste0(current.directory.scripts,'PNGs/Fig 5.png'),width=8, height=8,res=res,units=units)
#   }else{
#     pdf(paste0(current.directory.scripts,'PDFs/Fig 5.pdf'),width=8, height=8,paper='legal')
#   }
#   source('~/Desktop/Code/value.of.MSP.loop_interpol.R')
#   source('~/Desktop/Code/value.of.MSP.fx_2_interpol.R')
#   source('~/Desktop/Code/figure_5_code.R')
#   MSP.value.data=rbind(Mussel.value.tmp,Finfish.value.tmp,Kelp.value.tmp)
#   MSP.value.data.points=rbind(Mussel.value.tmp.points,Finfish.value.tmp.points,Kelp.value.tmp.points)
#   names(MSP.value.data)=c('Aquaculture.Value','Value.of.MSP','Group','Sector.Name','Sector.Type','Type.of.Conventional')
#   names(MSP.value.data.points)=c('Aquaculture.Value','Value.of.MSP','Group','Sector.Name','Sector.Type','Type.of.Conventional')
#   MSP.value.data$Sector.Name=factor(MSP.value.data$Sector.Name,levels=c('Mussel','Finfish','Kelp','Halibut','Viewshed','Benthic','Disease'))
#   Value.of.MSP.grid.plot<-ggplot(data=MSP.value.data)+
#     geom_point(data=subset(MSP.value.data[as.integer(MSP.value.data$Group)!=as.integer(MSP.value.data$Sector.Name),]),
#                aes(x=Aquaculture.Value,y=Value.of.MSP,shape=Type.of.Conventional,color=Type.of.Conventional),size=1.25)+
#     facet_grid(Sector.Name~Group)+scale_y_continuous(labels = percent,limits=c(0,1),breaks=c(0,.5,1))+
#     scale_x_continuous(labels = percent,breaks=c(0,.5,1))+
#     geom_text(data=subset(MSP.value.data[as.integer(MSP.value.data$Group)==as.integer(MSP.value.data$Sector.Name),]),x=.5,y=.5,size=15,label='NA')+
#     geom_point(data=MSP.value.data.points,aes(x=Aquaculture.Value,y=Value.of.MSP,shape=Type.of.Conventional,color=Type.of.Conventional),size=1.2)+
#     xlab("Aquaculture Value")+ylab("Value of Marine Spatial Planning")+
#     scale_colour_manual(name = "Conventional Planning :",labels=c("Constrained","Unconstrained"),values=c("coral1","mediumorchid1"))+
#     scale_shape_manual(name = "Conventional Planning :",labels=c("Constrained","Unconstrained"),values=c(16,24))+
#     scale_fill_manual(name = "Conventional Planning :",labels=c("Constrained","Unconstrained"),values=c("coral1","mediumorchid1"))+
#     #     geom_line(aes(x=c(0,1),y=c(0,0)),color="grey",linetype='dashed',size=1)+
#     theme_bw(base_size = 15)+theme(panel.grid=element_blank(),legend.position="bottom")
#   Value.of.MSP.grid.plot+theme(axis.title=element_text(size=12),panel.spacing = unit(1, "lines"),
#                                strip.background=element_rect(fill='white'),strip.text=element_text(size=10),axis.ticks.margin=unit(1,'lines'))
#   dev.off()
# }
#
# # ## Combine Figures
# # require(png)
# # require(grid)
# # png(paste0(current.directory.scripts,'/PNGs/MS_Figures.png'))
# # lapply(ll <- list.files(path = '~/Desktop/Code/MS Figures/PNGs',patt='.*[.]png'),function(x){
# #   img <- as.raster(readPNG(paste0('~/Desktop/Code/MS Figures/PNGs/',x)))
# #   grid.newpage()
# #   grid.raster(img, interpolate = FALSE)
#
# # })
# # dev.off()
#
#
#
#
# # plots <- lapply(ll <- list.files(path = paste0(current.directory.scripts,'PNGs/'), patt='.*[.]png'),function(x){
# #   img <- as.raster(readPNG(x))
# #   rasterGrob(img, interpolate = FALSE)
# # })
# # require(ggplot2)
# # require(gridExtra)
# # ggsave("multipage.pdf", marrangeGrob(grobs=plots, nrow=2, ncol=2))
#
#
#
# ## Figure 1
# ## SI Figures
# current.directory.code <- '~/Desktop/Code/SI Figures Code/'
# current.directory.figures <- '~/Desktop/Code/SI Figures/'
# # Coastline.data<-read.csv(file='Coastline.csv',header=F)
# theme_3 = theme(plot.margin = unit(c(.1,.1,.1,.1), units = "lines"),
#               axis.text = element_blank(),
#               axis.title = element_blank(),
#               axis.ticks = element_blank(),
#               # axis.ticks.length = unit(0, "lines"),
#               # axis.ticks.margin = unit(0, "lines"),
#               panel.background=element_rect(fill="white"),
#               panel.grid=element_blank(),
#               plot.title=element_text(hjust =.35,vjust=.1))
#
# theme_2 = theme(plot.margin = unit(c(.1,.1,.1,.1), units = "lines"),
#               axis.text = element_blank(),
#               axis.title = element_blank(),
#               axis.ticks = element_blank(),
#               # axis.ticks.length = unit(0, "lines"),
#               # axis.ticks.margin = unit(0, "lines"),
#               panel.background=element_rect(fill="white"),
#               panel.grid=element_blank(),
#               plot.title=element_text(hjust =.30,vjust=.1))
#
# labs = labs(x = NULL, y = NULL)
# # Figure S1
# for(itor in 1:2){
#   if(itor==1){
#     png(paste0(current.directory.figures,'PNGs/Fig S1.png'),width=width, height=height,res=res,units=units)
#   }else{
#     pdf(paste0(current.directory.figures,'PDFs/Fig S1.pdf'),width=width, height=height,paper='legal')
#   }
#   img_s1 <- readPNG(paste0(current.directory.code,'Fig S1.png'),native=T,info=T)
#   g_s1 <- rasterGrob(img_s1, interpolate=TRUE)
#
#   s1<-qplot(1:10, 1:10, geom="blank") +
#     annotation_custom(g_s1, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)+
#     theme_3 + labs
#   print(s1)
#   dev.off()}
# # Figure S2
# for(itor in 1:2){
#   if(itor==1){
#     png(paste0(current.directory.figures,'PNGs/Fig S2.png'),width=width, height=height,res=res,units=units)
#   }else{
#     pdf(paste0(current.directory.figures,'PDFs/Fig S2.pdf'),width=width, height=height,paper='legal')
#   }
#   img_S2A <- readPNG(paste0(current.directory.code,"S2A.png"),native=T,info=T)
#   g_S2A <- rasterGrob(img_S2A, interpolate=TRUE)
#
#   S2A<-qplot(1:10, 1:10, geom="blank") + ggtitle('A') +
#     annotation_custom(g_S2A, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)+
#     theme_3 + labs
#
#   img_S2B <- readPNG(paste0(current.directory.code,"S2B.png"),native=T,info=T)
#   g_S2B <- rasterGrob(img_S2B, interpolate=TRUE)
#
#   S2B<-qplot(1:10, 1:10, geom="blank") + ggtitle('B') +
#     annotation_custom(g_S2B, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)+
#     theme_3 + labs
#
#   img_S2C <- readPNG(paste0(current.directory.code,"S2C.png"),native=T,info=T)
#   g_S2C <- rasterGrob(img_S2C, interpolate=TRUE)
#
#   S2C<-qplot(1:10, 1:10, geom="blank") + ggtitle('C') +
#     annotation_custom(g_S2C, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)+
#     theme_3 + labs
#
#   grid.arrange(S2A,S2B,S2C,ncol=1,nrow=3,padding=unit(-.1,'line'))
#   dev.off()}
# # S3
#   for(itor in 1:2){
#     if(itor==1){
#       png(paste0(current.directory.figures,'PNGs/Fig S3.png'),width=width, height=height,res=res,units=units)
#     }else{
#       pdf(paste0(current.directory.figures,'PDFs/Fig S3.pdf'),width=width, height=height,paper='legal')
#     }
#   img_S3A <- readPNG(paste0(current.directory.code,"S3A.png"),native=T,info=T)
#   g_S3A <- rasterGrob(img_S3A, interpolate=TRUE)
#
#   S3A<-qplot(1:10, 1:10, geom="blank") + ggtitle('A') +
#     annotation_custom(g_S3A, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)+
#     theme_3 + labs
#
#   img_S3B <- readPNG(paste0(current.directory.code,"S3B.png"),native=T,info=T)
#   g_S3B <- rasterGrob(img_S3B, interpolate=TRUE)
#
#   S3B<-qplot(1:10, 1:10, geom="blank") + ggtitle('B') +
#     annotation_custom(g_S3B, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)+
#     theme_3 + labs
#
#   img_S3C <- readPNG(paste0(current.directory.code,"S3C.png"),native=T,info=T)
#   g_S3C <- rasterGrob(img_S3C, interpolate=TRUE)
#
#   S3C<-qplot(1:10, 1:10, geom="blank") + ggtitle('C') +
#     annotation_custom(g_S3C, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)+
#     theme_3 + labs
#   grid.arrange(S3A,S3B,S3C,ncol=1,nrow=3,padding=unit(-.1,'line'))
#   dev.off()}
# # S4
# for(itor in 1:2){
#   if(itor==1){
#     png(paste0(current.directory.figures,'PNGs/Fig S4.png'),width=width, height=height,res=res,units=units)
#   }else{
#     pdf(paste0(current.directory.figures,'PDFs/Fig S4.pdf'),width=width, height=height,paper='legal')
#   }
#   img_s4 <- readPNG(paste0(current.directory.code,"Fig S4.png"),native=T,info=T)
#   g_s4 <- rasterGrob(img_s4, interpolate=TRUE)
#
#   s4<-qplot(1:10, 1:10, geom="blank") +
#     annotation_custom(g_s4, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)+
#     theme_3 + labs
#   print(s4)
#   dev.off()}
# # S5
#   for(itor in 1:2){
#     if(itor==1){
#       png(paste0(current.directory.figures,'PNGs/Fig S5.png'),width=width, height=height,res=res,units=units)
#     }else{
#       pdf(paste0(current.directory.figures,'PDFs/Fig S5.pdf'),width=width, height=height,paper='legal')
#     }
#
#   img_S5A <- readPNG(paste0(current.directory.code,"S5A.png"),native=T,info=T)
#   g_S5A <- rasterGrob(img_S5A, interpolate=TRUE)
#
#   S5A<-qplot(1:10, 1:10, geom="blank") + ggtitle('A') +
#     annotation_custom(g_S5A, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)+
#     theme_2 + labs
#
#   img_S5B <- readPNG(paste0(current.directory.code,"S5B.png"),native=T,info=T)
#   g_S5B <- rasterGrob(img_S5B, interpolate=TRUE)
#
#   S5B<-qplot(1:10, 1:10, geom="blank") + ggtitle('B') +
#     annotation_custom(g_S5B, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)+
#     theme_2 + labs
#
#   grid.arrange(S5A,S5B,ncol=1,nrow=2,padding=unit(-.1,'line'))
#   dev.off()}
# # S6
# for(itor in 1:2){
#   if(itor==1){
#      png(paste0(current.directory.figures,'PNGs/Fig S6.png'),width=width, height=height,res=res,units=units)
#   }else{
#      pdf(paste0(current.directory.figures,'PDFs/Fig S6.pdf'),width=width, height=height,paper='legal')
#   }
#   img_S6A <- readPNG(paste0(current.directory.code,"S6A.png"),native=T,info=T)
#   g_S6A <- rasterGrob(img_S6A, interpolate=TRUE)
#
#   S6A<-qplot(1:10, 1:10, geom="blank") + ggtitle('A') +
#     annotation_custom(g_S6A, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)+
#     theme_2 + labs
#
#   img_S6B <- readPNG(paste0(current.directory.code,"S6B.png"),native=T,info=T)
#   g_S6B <- rasterGrob(img_S6B, interpolate=TRUE)
#
#   S6B<-qplot(1:10, 1:10, geom="blank") + ggtitle('B') +
#     annotation_custom(g_S6B, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)+
#     theme_2 + labs
#
#   grid.arrange(S6A,S6B,ncol=1,nrow=2,padding=unit(-.1,'line'))
#   dev.off()
# }
# # S7
# for(itor in 1:2){
#   if(itor==1){
#     png(paste0(current.directory.figures,'PNGs/Fig S7.png'),width=width, height=height,res=res,units=units)
#   }else{
#     pdf(paste0(current.directory.figures,'PDFs/Fig S7.pdf'),width=width, height=height,paper='legal')
#   }
#   img_S7A <- readPNG(paste0(current.directory.code,"S7A.png"),native=T,info=T)
#   g_S7A <- rasterGrob(img_S7A, interpolate=TRUE)
#
#   S7A<-qplot(1:10, 1:10, geom="blank") + ggtitle('A') +
#     annotation_custom(g_S7A, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)+
#     theme_3 + labs + theme(plot.title=element_text(hjust=0))
#
#   img_S7B <- readPNG(paste0(current.directory.code,"S7B.png"),native=T,info=T)
#   g_S7B <- rasterGrob(img_S7B, interpolate=TRUE)
#
#   S7B<-qplot(1:10, 1:10, geom="blank") + ggtitle('B') +
#     annotation_custom(g_S7B, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)+
#     theme_3 + labs + theme(plot.title=element_text(hjust=0))
#
#   img_S7C <- readPNG(paste0(current.directory.code,"S7C.png"),native=T,info=T)
#   g_S7C <- rasterGrob(img_S7C, interpolate=TRUE)
#
#   S7C<-qplot(1:10, 1:10, geom="blank") + ggtitle('C') +
#     annotation_custom(g_S7C, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)+
#     theme_3 + labs + theme(plot.title=element_text(hjust=0))
#
#   img_S7D <- readPNG(paste0(current.directory.code,"S7D.png"),native=T,info=T)
#   g_S7D <- rasterGrob(img_S7D, interpolate=TRUE)
#
#   S7D<-qplot(1:10, 1:10, geom="blank") + ggtitle('D') +
#     annotation_custom(g_S7D, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)+
#     theme_3 + labs + theme(plot.title=element_text(hjust=0))
#
#   grid.arrange(S7A,S7B,S7C,S7D,ncol=2,nrow=2,padding=unit(.5,'line'))
#   dev.off()
# }
# # S8
# for(itor in 1:2){
#   if(itor==1){
#     png(paste0(current.directory.figures,'PNGs/Fig S8.png'),width=width, height=height,res=res,units=units)
#   }else{
#     pdf(paste0(current.directory.figures,'PDFs/Fig S8.pdf'),width=width, height=height,paper='legal')
#   }
#   img_S8A <- readPNG(paste0(current.directory.code,"S8A.png"),native=T,info=T)
#   g_S8A <- rasterGrob(img_S8A, interpolate=TRUE)
#
#   S8A<-qplot(1:10, 1:10, geom="blank") + ggtitle('A') +
#     annotation_custom(g_S8A, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)+
#     theme_3 + labs
#
#   img_S8B <- readPNG(paste0(current.directory.code,"S8B.png"),native=T,info=T)
#   g_S8B <- rasterGrob(img_S8B, interpolate=TRUE)
#
#   S8B<-qplot(1:10, 1:10, geom="blank") + ggtitle('B') +
#     annotation_custom(g_S8B, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)+
#     theme_3 + labs
#
#   grid.arrange(S8A,S8B,ncol=1,nrow=2,padding=unit(.5,'line'))
#   dev.off()
# }
#
# # # setwd('C:/Users/Joel/Desktop/Thesis YTD/Code')
# #   # panel.EF<-function (x, y, itor=0, epsilon=.001, bg = NA, pch = 20, cex = .01, ...)
# #   # {
# #   #   x.MSP=x[color.vector=='coral']
# #   #   y.MSP=y[color.vector=='coral']
# #   #
# #   #   points(x.MSP,y.MSP, pch = 16, col = alpha("dodgerblue",1/75),cex = cex/2)
# #   #   x.U=x[color.vector=='purple']
# #   #   y.U=y[color.vector=='purple']
# #   #   x.S=x[color.vector=='green']
# #   #   y.S=y[color.vector=='green']
# #   #   x.EF=NULL
# #   #   y.EF=NULL
# #   #   alpha.mat.tmp=seq(from=0,by=epsilon,to=1)
# #   #   # MSP
# #   #   for(itor in 1:length(alpha.mat.tmp)){
# #   #     alpha.tmp=alpha.mat.tmp[itor]
# #   #     A=(alpha.tmp*x.MSP)+((1-alpha.tmp)*y.MSP)
# #   #     I=which(A==max(A))
# #   #     x.EF[itor]=max(unique(x.MSP[I]))
# #   #     I.tmp.x=which(x.MSP==max(unique(x.MSP[I])))
# #   #     I.tmp.y=which(y.MSP[I.tmp.x]==max(unique(y.MSP[I.tmp.x])))
# #   #     y.EF[itor]=unique(y.MSP[I.tmp.x[I.tmp.y]])}
# #   #   x.EF.original=x.EF;y.EF.original=y.EF;
# #   #   if(length(unique(x.EF.original))!=1&length(unique(x.EF.original))!=1){
# #   #     EF.inter=approx(x.EF.original,y=y.EF.original,n=length(alpha.mat.tmp))
# #   #     x.EF=EF.inter$x;y.EF=EF.inter$y;
# #   #   }else{
# #   #   }
# #   #   lines(sort(x.EF),y.EF[order(x.EF)],col="midnightblue",lwd=2,lty=1)
# #   #   lines(sort(x.U),y.U[order(x.U)],col = "mediumorchid1",lwd=2,lty=1)
# #   #   lines(sort(x.S),y.S[order(x.S)],col = "coral1",lwd=2,lty=1)
# #   # }
# #   #   pdf.options(width = 8, height = 6.4)
#
# #   # source('pairs2.R')
# # # for(itor in 1:2){
# # #   if(itor==1){
# #     png('Fig S9.png',width=8, height=6.4,res=res,units=units)
# #   # }else{
# #   #   pdf('Fig S9.pdf',width=8, height=6.4,paper='legal')
# #   # # }
# #   color.vector=color.vector.max # Color Vector For Seperating the MSP from Conventional Solutions
# #   pairs2(100*Master.matrix.max.pure.profit,lower.panel=panel.EF,
# #          upper.panel=NULL,col=color.vector,cex=0.8,xlim=c(0,100),
# #          ylim=c(0,100),pch=16,font.labels=3,cex.axis=1,las=1,xaxp=c(0,100,4),yaxp=c(0,100,2),
# #          gap=1)
=======
D.V <- rep(NA,times = length(F.V))
D.V[F.V > 0] <- readMat(paste0(wkdir,'/MSP_Model/Scripts/tmp.mat'))$d
system2('rm',args = paste0(wkdir,filename))
# Save all of the raw outputs of each sector model in a seperate file --> do later
print('Raw Impacts/Value.....')
Raw_Impacts <- data.frame(Mussel = M.V, Finfish = F.V, Kelp = K.V, Halibut = H.V,
  Viewshed_Mussel_Kelp = V_MK.V, Viewshed_Finfish = V_F.V, Benthic_Impacts = B.V,
  Disease_Risk = D.V) %>% glimpse()
# Make a .mat file of the sector files
writeMat('~/MSP_Model/Input/Data/Raw_Impacts.mat',Raw_Impacts = Raw_Impacts)
## Tradeoff Model
# Define parameters for the model
sector_names <- names(Raw_Impacts)
n <- 7 # Number of sectors
i <- 1061 # Number of sites, nrow(Raw_Impacts)
p <- 4 # Number of management options, 0 = no development, 1 = develop mussel, 2 = develop finfish, 3 = develop kelp, 4 = no development
# Using the definied parameters derive variable V_n_i_p (value to sector n at site i for pursuing development option p)
p_options <- list(c('Halibut'),c('Mussel','Viewshed_Mussel_Kelp'),c('Finfish','Viewshed_Finfish','Benthic_Impacts','Disease_Risk'),c('Kelp','Viewshed_Mussel_Kelp'))
# Make a list of dataframes consiting of the responses for each policy p for each sector n at each site i
R_n_i_p <- setNames(lapply(1:p,df = Raw_Impacts,FUN = function(x,df){
      # Make the dataframe for each policy
      tmp <- setNames(data.frame(t(do.call('rbind',lapply(1:length(sector_names), FUN = function(y){
        if(sector_names[y] %in% p_options[[x]]){
          return(df[[y]])
        }else{
          # For a sector recieving zero value or full impact set sector values to zero unless they were previously set as NA
          df[[y]][!is.na(df[[y]])] <- 0
          return(df[[y]])
        }
        })))),sector_names)
    }),c('No_Development','Develop_Mussel','Develop_Finfish','Develop_Kelp'))
# Add in the ifelse to make it comparable to lines 52 - 88 in CW code
R_n_i_p[[1]] <- R_n_i_p[[1]] %>%
  mutate(Viewshed = 0) %>%
  select(-Viewshed_Finfish,-Viewshed_Mussel_Kelp) %>% select(Mussel,Finfish,Kelp,Halibut,Viewshed,Benthic_Impacts,Disease_Risk)
R_n_i_p[[2]] <- R_n_i_p[[2]] %>%
  mutate(Viewshed = Viewshed_Mussel_Kelp) %>% mutate(Halibut = ifelse(Mussel == 0,Raw_Impacts$Halibut,0)) %>%
  mutate(Viewshed = ifelse(Mussel > 0,Raw_Impacts$Viewshed_Mussel_Kelp,0)) %>%
  select(-Viewshed_Finfish,-Viewshed_Mussel_Kelp) %>% select(Mussel,Finfish,Kelp,Halibut,Viewshed,Benthic_Impacts,Disease_Risk)
R_n_i_p[[3]] <- R_n_i_p[[3]] %>%
  mutate(Viewshed = Viewshed_Finfish) %>% mutate(Halibut = ifelse(Finfish == 0,Raw_Impacts$Halibut,0)) %>%
  mutate(Viewshed = ifelse(Finfish > 0,Raw_Impacts$Viewshed_Finfish,0)) %>%
  select(-Viewshed_Finfish,-Viewshed_Mussel_Kelp) %>% select(Mussel,Finfish,Kelp,Halibut,Viewshed,Benthic_Impacts,Disease_Risk)
R_n_i_p[[4]] <- R_n_i_p[[4]] %>%
  mutate(Viewshed = Viewshed_Mussel_Kelp) %>% mutate(Halibut = ifelse(Kelp == 0,Raw_Impacts$Halibut,0)) %>%
  mutate(Viewshed = ifelse(Kelp > 0,Raw_Impacts$Viewshed_Mussel_Kelp,0)) %>%
  select(-Viewshed_Finfish,-Viewshed_Mussel_Kelp) %>% select(Mussel,Finfish,Kelp,Halibut,Viewshed,Benthic_Impacts,Disease_Risk)
true_sector_names <- names(R_n_i_p[[1]])
# For sectors whose response is negative, calculate R_bar
R_negative_sector_names <- names(R_n_i_p[[1]])[grepl(c('Viewshed|Benthic_|Disease_'),names(R_n_i_p[[1]]))]
R_bar_n <- setNames(data.frame(t(do.call('rbind',lapply(1:length(R_negative_sector_names), data = lapply(R_n_i_p, R_negative_sector_names, FUN = function(x,y){
  apply(x[names(x) %in% R_negative_sector_names],MARGIN = 2, FUN = function(x) max(x,na.rm = T))
  }),FUN = function(x,data){
    max(sapply(data,"[",x),na.rm = T)
    })))),R_negative_sector_names)
# apply(sapply(1:1061,FUN = function(x){apply(data.frame(do.call("rbind",lapply(V_n_i_p,"[",x,))),MARGIN = 2,FUN = sum)}),MARGIN = 1,FUN = sum,na.rm=T)
# Then calculate V_n_i_p based on the response of each sector (Supp. Info, Eq. S26)
V_n_i_p <- setNames(lapply(1:p, df = R_n_i_p, df_bar = R_bar_n, FUN = function(x, df, df_bar){
  setNames(data.frame(t(do.call('rbind',lapply(1:length(names(R_n_i_p[[x]])), FUN = function(y){
    if(names(df[[x]])[y] %in% R_negative_sector_names){
      return(R_bar_n[[names(df[[x]])[y]]] - df[[x]][,y])
    }else{
      return(df[[x]][,y])
    }
    })))),names(R_n_i_p[[x]]))
  }), c('No_Development','Develop_Mussel','Develop_Finfish','Develop_Kelp'))
# apply(sapply(1:7,FUN = function(x){apply(data.frame(do.call("rbind",lapply(V_n_i_p,"[",,x))),MARGIN = 2,FUN = max)}),MARGIN = 1,FUN = sum,na.rm=T)
# Scale sector value n for each site i, by the domain-wide value that would be attained if the ideal development option to the sector was selected at site i, X_n_i_p (Supp. Info, Eq. S27)
# X_n_i_p <- setNames(lapply(1:p, df = V_n_i_p, FUN = function(p,df){
#       return(setNames(data.frame(t(do.call('rbind',lapply(1:length(sector_names), FUN = function(n){
#         if((p == 2 | p == 4) & sector_names[n] == 'Viewshed_Finfish'){
#           return(rep(0,times = i))
#         }else if((p == 3 | p == 1) & sector_names[n] == 'Viewshed_Mussel_Kelp'){
#           return(rep(0,times = i))
#         }else{
#           if(p == 1 & sector_names[n] == 'Viewshed_Finfish'){
#             return(apply(df[[p]][,c(5,6)],MARGIN = 1, FUN = max,na.rm = T) / sum(apply(data.frame(apply(df[[p]][,c(5,6)],MARGIN = 1, FUN = max,na.rm = T)),MARGIN = 1, FUN = function(z){ifelse(!all(is.na(z)),max(z, na.rm = T),NA)}),na.rm = T))
#           }else{
#             return(df[[p]][,n] / sum(apply(sapply(df,"[", ,n),MARGIN = 1, FUN = function(z){ifelse(!all(is.na(z)),max(z, na.rm = T),NA)}),na.rm = T))
#           }
#         }
#       })))),true_sector_names))
#   }),c('No_Development','Develop_Mussel','Develop_Finfish','Develop_Kelp'))
# sum(apply(sapply(V_n_i_p,"[", ,7),MARGIN = 1, FUN = function(z){ifelse(!all(is.na(z)),max(z, na.rm = T),NA)}),na.rm = T)
# sum(apply(sapply(V_n_i_p,"[", ,6),MARGIN = 1, FUN = function(z){ifelse(!all(is.na(z)),max(z, na.rm = T),NA)}),na.rm = T)
X_n_i_p <- setNames(lapply(1:p, df = V_n_i_p, FUN = function(p,df){
      return(setNames(data.frame(t(do.call('rbind',lapply(1:length(true_sector_names), FUN = function(n){
        print(n);
        print(p);
        print(sum(apply(sapply(df,"[", ,n),MARGIN = 1, FUN = function(z){ifelse(!all(is.na(z)),max(z, na.rm = T),NA)}),na.rm = T));
        df[[p]][,n] / sum(apply(sapply(df,"[", ,n),MARGIN = 1, FUN = function(z){ifelse(!all(is.na(z)),max(z, na.rm = T),NA)}),na.rm = T)
      })))),true_sector_names))
  }),c('No_Development','Develop_Mussel','Develop_Finfish','Develop_Kelp'))
# Create the sector weights (alpha's)
library(gtools)
epsilon <- .20 # Epsilon step size, default is 0.20
a_values <- seq(from = 0, to = 1, by = epsilon) # The unique values for each sector and site
# a_tmp <- permutations(n = length(a_values),7,a_values,repeats.allowed=T) # Matrix of unique alpha weights
# a <- cbind(a_tmp[,1:5],a_tmp[,5],a_tmp[,6:7]) # Do not do this if truely only 7 sectors are used, i.e. viewshed maximum across both sets is used vs. how it is currently where they are handled seperately
a <- permutations(n = length(a_values),7,a_values,repeats.allowed=T)
# Find the optimal policy option for each site in a given, alpha
if(readline('Perform Full Analysis? Y/N ') == 'Y'){
  print('Finding optimal solutions.................')
  print_a <- seq(from = 0, to = nrow(a), by = 10000)
  obj_i <- sapply(1:nrow(a), FUN = function(x){
    if(x %in% print_a){print(paste0(x,' iterations'))}
    apply(sapply(1:p, df = X_n_i_p, FUN = function(y,df){
      apply(data.frame(mapply('*',df[[y]],c(a[x,]),SIMPLIFY = FALSE)), MARGIN = 1, FUN = function(z) sum(z,na.rm = T)) # Multiply each i for a given p by the sector specific weight set by a given row of the alpha matrix
      }),MARGIN = 1, which.max) - 1
  })
  x = 55988
  y = 1
  i.test <- 1000
  sapply(55988, FUN = function(x){
    if(x %in% print_a){print(paste0(x,' iterations'))}
    apply(sapply(1:p, df = lapply(X_n_i_p, "[",i.test,), FUN = function(y,df){
      apply(data.frame(mapply('*',df[[y]],c(a[x,]),SIMPLIFY = FALSE)), MARGIN = 1, FUN = function(z) sum(z,na.rm = T)) # Multiply each i for a given p by the sector specific weight set by a given row of the alpha matrix
      }),MARGIN = 1, which.max) - 1
  })
  # # Save model results
  write.table(x = data.frame(obj_i,stringsAsFactors = F),file = file.path(paste0(wkdir,'/MSP_Model/Output/Data/MSP_Planning_Results.csv')), sep = ",",quote = FALSE, col.names = FALSE, row.names = FALSE)
  ## Turn primary variables into a list
  JS_MSP.list <- setNames(list(Raw_Impacts,R_n_i_p,R_bar_n,V_n_i_p,X_n_i_p,a),c('Raw_Impacts','R_n_i_p','R_bar_n','V_n_i_p','X_n_i_p','a'))
  # Save primary variables to load into R
  save(JS_MSP.list,file = '~/MSP_Model/Output/Data/JS_MSP_Output.Rdata')
  # Save primary variables to load into Matlab
  writeMat('~/MSP_Model/Output/Data/JS_MSP_Output.mat',Raw_Impacts = Raw_Impacts, R_n_i_p = R_n_i_p, R_bar_n = R_bar_n, V_n_i_p = V_n_i_p, X_n_i_p = X_n_i_p, a = a)
}else{
  print('loading planning results')
  obj_i <- read.csv(file.path('~/MSP_Model/Output/Data/MSP_Planning_Results.csv'))
}
# Summarize JS Results
JS_Summary <- do.call('rbind',lapply(1:i, FUN = function(x){
                print(x)
                # return(table(obj_i[x,]))
                tmp <- rep(0,times = p)
                summary.tmp <- table(obj_i[x,])
                choices <- as.numeric(names(summary.tmp))
                numbers <- unname(summary.tmp)
                tmp[choices + 1] <- numbers
                return(tmp)
              }))
# Load CW Data
CW_Summary <- read.csv('~/MSP_Model/Scripts/CrowT0v1/CW_Results_Summary.csv',header = FALSE) %>% glimpse()
# CW_Summary <- readMat('~/MSP_Model/Scripts/CrowT0v1/CW_Results_Summary.mat')[['CW.results']]
for(itor in 1:100){
  test.i <- sample(1:1061,1)
  print(test.i)
  print(CW_Summary[test.i,])
  print(JS_Summary[test.i,])
  print(CW_Summary[test.i,] - JS_Summary[test.i,])
}


summary(JS_Summary - CW_Summary)
setNames(apply(JS_Summary - CW_Summary,MARGIN = 2, FUN = function(x){length(which(x != 0))}),c('No_Development','Develop_Mussel','Develop_Finfish','Develop_Kelp'))
foo <- setNames(sapply(1:p,FUN = function(x){
  which(JS_Summary[[x]] - CW_Summary[[x]] != 0)
  }),c('No_Development','Develop_Mussel','Develop_Finfish','Develop_Kelp'))
sapply(1:i,FUN = function(x){JS_Summary[x,] - CW_Summary[x,] != 0})
test.i <- 1000
print('Raw_Impacts')
print(Raw_Impacts[test.i,])
print('................................................................................................')
print('R_n_i_p')
print(sapply(R_n_i_p,"[",test.i,))
print('................................................................................................')
print('R_bar_n')
print(R_bar_n)
print('................................................................................................')
print('V_n_i_p')
print(sapply(V_n_i_p,"[",test.i,))
print('................................................................................................')
print('X_n_i_p')
print(sapply(X_n_i_p,"[",test.i,))
CW_Variables <- readMat('~/MSP_Model/Scripts/CrowT0v1/TOA_data.mat')
X_n_i_p.CW = setNames(lapply(1:p,FUN = function(x){setNames(data.frame(CW_Variables[['X.n.i.p']][,,x]),true_sector_names)}),c('No_Development','Develop_Mussel','Develop_Finfish','Develop_Kelp'))
lapply(1:p,FUN = function(x){
  X_n_i_p.CW[[x]] - X_n_i_p[[x]]
  })
X_n_i_p.diff <- setNames(lapply(1:p,FUN = function(x){setNames(data.frame(round(X_n_i_p.CW[[x]] - X_n_i_p[[x]],10)),true_sector_names)}),c('No_Development','Develop_Mussel','Develop_Finfish','Develop_Kelp'))
apply(rbindlist(X_n_i_p.diff),MARGIN = 2,FUN = sum,na.rm = T)
# source('~/MSP_Model/Scripts/Model_QA.r')
# Load and subtract 1 from each entry in order to allow for an accurate comparison
# obj_i.JS <- obj_i
# obj_i.CW <- apply(read.csv('~/MSP_Model/Scripts/CrowT0v1/Policy_i_a.csv',stringsAsFactors=FALSE, header=T, nrows=1061),MARGIN = 2,FUN = function(x){x - 1})
# # Compare the two sets of results
# obj_i.compare <- obj_i.JS == obj_i.CW
# # Summarize for each row
# summary.i <- apply(obj_i.compare, MARGIN = 1, FUN = function(x){length(which(x))})
# summary.j <- apply(obj_i.compare, MARGIN = 2, FUN = function(x){length(which(x))})
#
# # Compare the two case study
# i.test <- 138
# a.test <- 55988
#
# print(obj_i.JS[i.test,a.test] == obj_i.CW[i.test,a.test])


# # Add the actual values of each sector n for each policy p at site i
# # Develop mussel
# p = p_options[[1]]
# # Responses
# # No development
# i = 1
# p = 0
# n_v <- 'Halibut'
# n_0 <- sector_names[!(sector_names %in% names_v)]
# str(lapply(1:p,df = Raw_Impacts,FUN = Responses))
# df = Raw_Impacts
# Responses <- function(x,df){
#   return(t(do.call('rbind',lapply(1:length(sector_names), FUN = function(y){
#     if(sector_names[y] %in% p_options[[x]]){
#       return(df[[y]])
#     }else{
#       df[[y]][!is.na(df[[y]])] <- 0
#       return(df[[y]])
#     }
#     }))))
# }
#   df[,p_options[[x]]]
#   tdf[,!(sector_names %in% p_options[[x]])]
#   unlist(ifelse(sector_names %in% p_options[[x]],df %>% select(contains(p_options[[x]])),rep(0,times = nrow(df))))
#   df[,!(sector_names %in% p_options[[x]])] <- ifelse(sector_names %in% p_options[[x]])
#
#
# }
# sector_names[!(sector_names %in% p_options[[x]])]
#
# # First invert the impacts using Eq. 26 of the Supplemental Information in Lester et al. 2017
# Response <- function(R,n){
#
#   if(Sector.Name %in% c('Mussel','Finfish','Kelp','Halibut')){
#     Response <- Value
#   }else{
#     if(Impacted == T){
#       Inverted <- max(Value) - Value
#     }else{
#       Inverted = Value
#     }
#     # Scale each inverted by the max value
#     Inverted.Scaled <- Inverted / max(Inverted)
#     Non.Inverted.Scaled <- Value / sum(Value)
#   }
#   # Logical for places which have zero value
#   # Logical.Value <- 1 * (Value>0)
#   # Make a data frame for value, inverted value, and inverted scaled value
#   X.df <- data.frame(Value = Value, Inverted = Inverted, Inverted_Scaled = Inverted.Scaled, Non_Inverted_Scaled = Non.Inverted.Scaled)
#   return(X.df)
# }
#
# Inverted <- rep(NA,length(Value))
# Inverted.Scaled <- rep(NA,length(Value))
# Non.Inverted.Scaled <- rep(NA,length(Value))
# Value.tmp <- rep(NA,length(Value))
#
# Value <- Value[F.NPV[Aqua.Full.Domain.Logical]>0]
# Inverted[F.NPV[Aqua.Full.Domain.Logical]>0] <- max(Value) - Value
#
# Inverted.Scaled <- Inverted / max(Inverted,na.rm = T)
# Non.Inverted.Scaled[F.NPV[Aqua.Full.Domain.Logical]>0] <- Value / sum(Value,na.rm=T)
#
# Value.tmp[F.NPV[Aqua.Full.Domain.Logical]>0] <- Value
# Value <- Value.tmp
#
#
#
#
#
#
# Vi <- apply(cbind(F.Viewshed,MK.Viewshed),MARGIN = 1, FUN = max)[Aqua.Full.Domain.Logical]
# Vi_diff <- Vi - MK.Viewshed[Aqua.Full.Domain.Logical]
#
# New.Approach.DF <- function(Value,Impacted,Sector.Name){
#   # max subtracted from each value
#   if(Sector.Name == 'Benthic' | Sector.Name == 'Disease'){
#     Inverted <- rep(NA,length(Value))
#     Inverted.Scaled <- rep(NA,length(Value))
#     Non.Inverted.Scaled <- rep(NA,length(Value))
#     Value.tmp <- rep(NA,length(Value))
#
#     Value <- Value[F.NPV[Aqua.Full.Domain.Logical]>0]
#     Inverted[F.NPV[Aqua.Full.Domain.Logical]>0] <- max(Value) - Value
#
#     Inverted.Scaled <- Inverted / max(Inverted,na.rm = T)
#     Non.Inverted.Scaled[F.NPV[Aqua.Full.Domain.Logical]>0] <- Value / sum(Value,na.rm=T)
#
#     Value.tmp[F.NPV[Aqua.Full.Domain.Logical]>0] <- Value
#     Value <- Value.tmp
#   }else{
#     if(Impacted == T){
#       Inverted <- max(Value) - Value
#     }else{
#       Inverted = Value
#     }
#     # Scale each inverted by the max value
#     Inverted.Scaled <- Inverted / max(Inverted)
#     Non.Inverted.Scaled <- Value / sum(Value)
#   }
#   # Logical for places which have zero value
#   # Logical.Value <- 1 * (Value>0)
#   # Make a data frame for value, inverted value, and inverted scaled value
#   X.df <- data.frame(Value = Value, Inverted = Inverted, Inverted_Scaled = Inverted.Scaled, Non_Inverted_Scaled = Non.Inverted.Scaled)
#   return(X.df)
# }
# V.M.DF <- New.Approach.DF(Aqua$M.NPV,F,'Mussel')
# V.F.DF <- New.Approach.DF(Aqua$F.NPV,F,'Finfish')
# V.K.DF <- New.Approach.DF(Aqua$K.NPV,F,'Kelp')
# V.H.DF <- New.Approach.DF(H.NPV,F,'Halibut')
# V.V.DF <- New.Approach.DF(Vi,T,'Viewshed')
# V.V.Diff.DF <- New.Approach.DF(Vi_diff,T,'Viewshed')
# V.B.DF <- New.Approach.DF(Bi,T,'Benthic')
# V.D.DF <- New.Approach.DF(Di,T,'Disease')
# head(V.D.DF)
# # Combine Data into list
# library(gtools)
# a_values <- seq(from = 0, to = 1, by = epsilon)
# a <- permutations(n = length(a_values),7,a_values,repeats.allowed=T)
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
# # Calculate R_max for viewshed for each developable site
# Vi <- apply(cbind(F.Viewshed,MK.Viewshed),MARGIN = 1, FUN = max)[Aqua.Full.Domain.Logical]
# # Find the maximum response across all sites
# V.R_max <- max(Vi)
# # Invert Finfish responses and Mussel/Kelp responses
# VF.V_n_i_p <- V.R_max - F.Viewshed[Aqua.Full.Domain.Logical]
# VMK.V_n_i_p <- V.R_max - MK.Viewshed[Aqua.Full.Domain.Logical]
# ## Load Benthic Data
# Bi <- df$TOC.flux[F.NPV > 0]
# # Find the maximum response across all sites
# B.R_max <- max(Bi)
# # Invert Finfish responses and Mussel/Kelp responses
# B.V_n_i_p <- rep(0, times = length(which(Aqua.Full.Domain.Logical)))
# B.V_n_i_p[F.NPV[Aqua.Full.Domain.Logical] > 0] <- B.R_max - Bi
# # Disease
# ## Load disease connectivity matrix, which will be used as the adjacency matrix
# # for the disease propagation network
# disease_mat <- readMat(paste0(wkdir,'/MSP_Model/Input/Data/disease_connect_matrix.mat'))$disease.connect.matrix
# # Remove all sites which cannot be developed for finfish
# disease_mat_finfish <- disease_mat[F$Annuity > 0,F$Annuity > 0]
# system2('/Applications/MATLAB_R2016b.app/bin/matlab',
#   args = c('-nodesktop','-noFigureWindows','-nodisplay','-nosplash',
#   '-r \\"try, r'))
# graph <- graph_from_adjacency_matrix(disease_mat_finfish, weighted = T, mode = 'undirected')
# Di_compare <- eigen_centrality(graph, directed = F, weights = E(graph)$weight)$vector
# # plot(1:392,Di_compare/sum(Di_compare))
# # points(1:392,Di/sum(Di),col='red')
# Di <- read.csv(file = paste0(model_directory,'Data/Raw_Patch_Data.csv'))[F.NPV[Aqua.Full.Domain.Logical] > 0,7]
# # Find the maximum response across all sites
# D.R_max <- max(Di)
# # Invert Finfish responses and Mussel/Kelp responses
# D.V_n_i_p <- rep(0, times = length(which(Aqua.Full.Domain.Logical)))
# D.V_n_i_p[F.NPV[Aqua.Full.Domain.Logical] > 0] <- D.R_max - Di
# ## Load Benthic Data
# ## Disease Model
#
# ## Invert impacted sector values
# ## Scale annuities for all sectors
# ## Run MSP Model
# ## Run Dynamic Halibut Model
# ## Plot Results
#
#
#
#
#
#
#
#
#
# # run_matlab_script(paste0(wkdir,'/MSP_Model/Scripts/Halibut_Model_2016/Halibut_tuner_free_params_v4.m'))
# #
# # # Insert the directory in which the model folder is located
# # model_directory <- '~/Desktop/Aquaculture_MSP_Model/'
# # # Insert directory in which MATLAB is located
# # matlab_directory <- '/Applications/MATLAB_R2016b.app/bin'
# #
# # # Run MSP tradeoff model in MATLAB
# # run_matlab_script('~/Desktop/Aquaculture_Paper/Code/TOA_AquaMSP_CrowCode_v1NaN.m')
# # run_matlab_script('~/Desktop/Aquaculture_Paper/Code/EFpayoffs_AquaMSP_CrowCode_v2.m')
# #
# # n.sector <- 7 # Number of sectors
# # epsilon <- 0.2 # Stepsize of sector weights
# #
# # t <- 10 # Time Horizon
# # r <- 0.05 # Discount rate
# #
# # # Load sector data calculated from sector-specific bioeconomic models
# #   df <- read.csv(file=paste0(model_directory,'Data/SeaGrant_data_complete_2015.csv'),stringsAsFactors=FALSE)
# #   fulldomain <- df$TARGET_FID # Model domain
# #   discount_factor <- 1/((1+r)^c(1:t))
# #   r_iy_aqua <- do.call(rbind, replicate(length(fulldomain), discount_factor, simplify=FALSE))
# # #### Calculate Responses (R) to Development Options
# # ## Sectors which higher responses increase value
# #   # Aquaculture
# #     Value.KM <- function(yield,upfront.cost,annual.cost,price){
# #       revenue <- do.call(cbind, replicate(t, yield * price, simplify=FALSE))
# #       cost <- cbind(upfront.cost + annual.cost,
# #         do.call(cbind, replicate(t - 1, annual.cost, simplify=FALSE)))
# #       profit <- (revenue - cost) * r_iy_aqua
# #       profit[profit < 0] <- 0
# #       NPV <- apply(profit, FUN = sum, MARGIN = 1)
# #       Annuity = (r*NPV)/(1-((1+r)^-t))
# #       return(list(NPV = NPV,Annuity = Annuity))
# #     }
# #     Value.F <- function(yield,costs,price){
# #       revenue <- do.call(cbind, replicate(t, yield * price, simplify=FALSE))
# #       cost <- df$fish.annual.operating.costs
# #       profit <- (revenue - cost) * r_iy_aqua
# #       profit[profit < 0] <- 0
# #       NPV <- apply(profit, FUN = sum, MARGIN = 1)
# #       Annuity = (r*NPV)/(1-((1+r)^-t))
# #       return(list(NPV = NPV,Annuity = Annuity))
# #     }
# #
# #     M <- Value.KM(df$mussel.yield,
# #       df$mussel.upfront.cost,
# #       df$mussel.annual.operating.cost,3.3)
# #     F <- Value.F(df$fish.yield,
# #       df$fish.upfront.cost,
# #       unique(df$fish.price[df$fish.price>0]))
# #     K <- Value.KM(df$kelp.yield,
# #       df$kelp.upfront.cost,
# #       df$kelp.annual.operating.cost,3)
# #
# #     # Remove unprofitable cells to give the raw
# #     Aqua.Full.Domain <- data.frame(M$Annuity,F$Annuity,K$Annuity)
# #     Aqua.Full.Domain.Logical <- apply(1 * (Aqua.Full.Domain > 0),FUN=sum,MARGIN=1) > 0
# #     # Aqua <- Aqua.Full.Domain[Aqua.Full.Domain.Logical,]
# #     M.V_n_i_p <- M$Annuity[Aqua.Full.Domain.Logical]
# #     F.V_n_i_p <- F$Annuity[Aqua.Full.Domain.Logical]
# #     K.V_n_i_p <- K$Annuity[Aqua.Full.Domain.Logical]
# #   # Halibut
# #   Halibut_10yrNPVi <- read.csv(file=paste0(model_directory,'Data/Target_FID_and_Yi_fulldomain_NPV_at_MSY_noAqua.csv'));
# #   names(Halibut_10yrNPVi) <- c('FID','Yi')
# #   H.NPV <- Halibut_10yrNPVi[Aqua.Full.Domain.Logical,2]
# #   H.V_n_i_p <- H.NPV
# # ## Sectors which higher responses decreases value
# #   # Viewshed
# #     F.Viewshed <- as.numeric(gsub(",", "", df$res_views_8k)) + as.numeric(gsub(",", "",df$park_views_8k))
# #     MK.Viewshed <- as.numeric(gsub(",", "", df$res_views_3k)) + as.numeric(gsub(",", "",df$park_view_3k))
# #     # Calculate R_max for viewshed for each developable site
# #     Vi <- apply(cbind(F.Viewshed,MK.Viewshed),MARGIN = 1, FUN = max)[Aqua.Full.Domain.Logical]
# #     # Find the maximum response across all sites
# #     V.R_max <- max(Vi)
# #     # Invert Finfish responses and Mussel/Kelp responses
# #     VF.V_n_i_p <- V.R_max - F.Viewshed[Aqua.Full.Domain.Logical]
# #     VMK.V_n_i_p <- V.R_max - MK.Viewshed[Aqua.Full.Domain.Logical]
# #   # Benthic
# #     Bi <- df$TOC.flux[F.NPV > 0]
# #     # Find the maximum response across all sites
# #     B.R_max <- max(Bi)
# #     # Invert Finfish responses and Mussel/Kelp responses
# #     B.V_n_i_p <- rep(0, times = length(which(Aqua.Full.Domain.Logical)))
# #     B.V_n_i_p[F.NPV[Aqua.Full.Domain.Logical] > 0] <- B.R_max - Bi
# #   # Disease
# #     ## Load disease connectivity matrix, which will be used as the adjacency matrix
# #     # for the disease propagation network
# #     disease_mat <- readMat(paste0(model_directory,'Data/disease_connect_matrix.mat'))$disease.connect.matrix
# #     # Remove all sites which cannot be developed for finfish
# #
# #     ## 12/29/16 differences between Matlab and R eigenvector centrality metrics
# #     disease_mat_finfish <- disease_mat[F$Annuity > 0,F$Annuity > 0]
# #     graph <- graph_from_adjacency_matrix(disease_mat_finfish,diag = T, weighted = T, mode = 'undirected')
# #     Di_compare <- eigen_centrality(graph, directed = T, weights = E(graph)$weight)$vector
# #     # plot(1:392,Di_compare/sum(Di_compare))
# #     # points(1:392,Di/sum(Di),col='red')
# #     Di <- read.csv(file = paste0(model_directory,'Data/Raw_Patch_Data.csv'))[F.NPV[Aqua.Full.Domain.Logical] > 0,7]
# #     # Find the maximum response across all sites
# #     D.R_max <- max(Di)
# #     # Invert Finfish responses and Mussel/Kelp responses
# #     D.V_n_i_p <- rep(0, times = length(which(Aqua.Full.Domain.Logical)))
# #     D.V_n_i_p[F.NPV[Aqua.Full.Domain.Logical] > 0] <- D.R_max - Di
# # ### Scale Sectors
# #   Scaling <- function(Sector){Scaled_Sector <- Sector / max(Sector)}
# #    M.X <- Scaling(M.V_n_i_p)
# #    F.X <- Scaling(F.V_n_i_p)
# #    K.X <- Scaling(K.V_n_i_p)
# #    H.X <- Scaling(H.V_n_i_p)
# #    VF.X <- Scaling(VF.V_n_i_p)
# #    VMK.X <- Scaling(VMK.V_n_i_p)
# #    D.X <- Scaling(D.V_n_i_p)
# #
# #   Scaled_Data <- data.frame(M_X = M.X, F_X = F.X, K_X = K.X,
# #     H_X = H.X, VF_X = VF.X, VMK_X = VMK.X, D_X = D.X)
# # ## Analyze data to make extract all values
# # # Load old version results
# # raw_values.df <- read.csv('~/Desktop/CrowTOv1/Raw_Patch_Data.csv')
# # names(raw_values.df) <- c('M','F','K','H','V_F','V_MK','B','D')
# # write.csv(file = '~/Desktop/Aquaculture_Paper/Code/Raw_Patch_Data.csv',raw.annuities)
# #
# # data.df <- read.csv(file='~/Desktop/Code/MSP Planning Results April 2016/Dynamic_Values_Export.csv',header=F)
# # names(data.df) <- c('M','F','K','H','V','B','D')
# # plans.df <- read.csv(file='~/Desktop/Code/MSP Planning Results April 2016/Static_plans.csv',header=F)
# #
# # # Conventional Models
# # Unconstrained_data <- read.csv(file='~/Desktop/Code/MSP Planning Results April 2016/Unconstrained_Dynamic_Values_April.csv',header=F)
# # names(Unconstrained_data)=c('Mussel','Finfish','Kelp','Halibut','Viewshed','Benthic','Disease')
# # Constrained_data <- read.csv(file='~/Desktop/Code/MSP Planning Results April 2016/Constrained_Dynamic_Values_April.csv',header=F)
# # names(Constrained_data)=c('Mussel','Finfish','Kelp','Halibut','Viewshed','Benthic','Disease')
# #
# # monetry_values.df <- data.frame(M = data.df$M * sector_totals.df[1],
# #   F = data.df$F * sector_totals.df[2],
# #   K = data.df$K * sector_totals.df[3],
# #   H = data.df$H * sector_totals.df[4])
# # # Conversion of NPV to equivalent annuity
# # # Calculate annuities
# # annuities <- monetry_values.df %>% mutate(M = (r*M)/(1-((1+r)^-t)),
# #                 F = (r*F)/(1-((1+r)^-t)),
# #                   K = (r*K)/(1-((1+r)^-t)),
# #                   H = (r*H)/(1-((1+r)^-t)))
# # # Calculate the annuitiy for > 25% of Mussel and loss of annuity of less than 1%
# # data.df %>% filter(M > .25, 1 - H < .01) %>% select(M,H) %>%
# #       mutate(M = M * sector_totals.df[1],H = (1-H) * sector_totals.df[4])
# #
# # # Calculate the case study
# # Case_Study.df <- data.df %>% mutate(ID = 1:nrow(data.df)) %>% filter(M >= .05, F >= .05, K >= .05, H >= .95, V >= .95, B >= .95, D >= .95)
# # Case_Study.df[576,] %>% select(M,F,K) %>% mutate(M = M * sector_totals.df[1], F = F * sector_totals.df[2], K = K * sector_totals.df[3])
# # summary(factor(plans.df[,Case_Study.df$ID[576]]))
# #
# # # Calculate sector differences
# # color.vector=c(rep('coral',length.out=nrow(data.df)),rep('purple',length.out=1061),rep('green',length.out=1061))
# #  panel.EF<-function(x, y, itor=0, epsilon=.001, bg = NA, pch = 20, cex = .01, ...)
# #   {
# #     x.MSP=x[color.vector=='coral']
# #     y.MSP=y[color.vector=='coral']
# #
# #     # points(x.MSP,y.MSP, pch = 16, col = alpha("lightblue1",1/100),cex = cex)
# #     x.U=x[color.vector=='purple']
# #     y.U=y[color.vector=='purple']
# #     x.S=x[color.vector=='green']
# #     y.S=y[color.vector=='green']
# #     x.EF=NULL
# #     y.EF=NULL
# #     alpha.mat.tmp=seq(from=0,by=epsilon,to=1)
# #     # MSP
# #     for(itor in 1:length(alpha.mat.tmp)){
# #       alpha.tmp=alpha.mat.tmp[itor]
# #       A=(alpha.tmp*x.MSP)+((1-alpha.tmp)*y.MSP)
# #       I=which(A==max(A))
# #       x.EF[itor]=max(unique(x.MSP[I]))
# #       I.tmp.x=which(x.MSP==max(unique(x.MSP[I])))
# #       I.tmp.y=which(y.MSP[I.tmp.x]==max(unique(y.MSP[I.tmp.x])))
# #       y.EF[itor]=unique(y.MSP[I.tmp.x[I.tmp.y]])}
# #     x.EF.original=x.EF;y.EF.original=y.EF;
# #     if(length(unique(x.EF.original))!=1&length(unique(x.EF.original))!=1){
# #       EF.inter=approx(x.EF.original,y=y.EF.original,n=length(alpha.mat.tmp))
# #       x.EF=EF.inter$x;y.EF=EF.inter$y;
# #     }else{
# #     }
# #     # lines(sort(x.EF),y.EF[order(x.EF)],col="midnightblue",lwd=2,lty=1)
# #     # lines(sort(x.U),y.U[order(x.U)],col = "mediumorchid1",lwd=2,lty=1)
# #     # lines(sort(x.S),y.S[order(x.S)],col = "coral1",lwd=2,lty=1)
# #     return(cbind(x.EF,y.EF))
# #   }
# #
# # M_H.EF <- data.frame(panel.EF(data.df$M,data.df$H)) %>%
# #     mutate(M_raw.EF = x.EF * sector_totals.df[1],
# #     H_raw.EF = y.EF * sector_totals.df[4])
# # M_H.S <- data.frame(Constrained_data[,1], Constrained_data[,4], Constrained_data[,1] * sector_totals.df[1], Constrained_data[,4] * sector_totals.df[4])
# # names(M_H.S) <- c('M.percentage', 'H.percentage', 'M.annuity', 'H.annuity')
# #
# # M_H.EF %>% filter(x.EF == .440)
# # tail(M_H.S)
# #
# #
# # # Load Scenario Data
# # # print('Loading Data........')
# # # V1_JS.data <- read.csv(file = '~/Desktop/Code/V1_Policy_JS.csv',header = F)
# # # V2_JS.data <- read.csv(file = '~/Desktop/Code/V2_Policy_JS.csv',header = F)
# # # V2_CW.data <- read.csv(file = '~/Desktop/Code/V2_Policy_CW.csv',header = F)
# #
# # # freq_scenario <- function(data){
# # #   List.Data <- list()
# # #   fx.t <- proc.t()
# # #   print('This may take a while...........')
# # #   for(index in 1:ncol(data)){
# # #     if(index %in% seq(from = 0, to = ncol(data), by = 10000)){
# # #       print(index)
# # #       ptm <- proc.t()
# # #     }
# # #     tab.tmp <- data.frame(table(factor(data[,index],
# # #       levels = c(1:4),
# # #       labels = c('No Development','Mussel','Finfish','Kelp'))))
# # #     names(tab.tmp) <- c('Policy','Frequency')
# # #     list.name <- paste("Scenario_", index, sep="")
# # #     List.Data[[list.name]] <- tab.tmp %>% mutate(Scenario = rep(index,ts = length(tab.tmp)))
# # #   }
# # #   new.data <- rbindlist(List.Data)
# # #   if(index %in% seq(from = 0, to = ncol(data), by = 10000)){
# # #      return(print(proc.t()-ptm))
# # #   }
# # #   print(paste('t elapsed',proc.t() - fx.t))
# # #   return(new.data)
# # # }
# # # print('Arranging Data........')
# # # SW_Version.df <- freq_scenario(V1_JS.data) %>% mutate(Version = 'Old Approach')
# # # JS_Version.df <- freq_scenario(V2_JS.data) %>% mutate(Version = 'New Approach JS')
# # # CW_Version.df <- freq_scenario(V2_CW.data) %>% mutate(Version = 'New Approach CW')
# #
# # # # Combine all
# # # Version.df <- bind_rows(SW_Version.df %>% mutate(Version = 'Old Approach'),
# # #     JS_Version.df %>% mutate(Version = 'New Approach JS'),
# # #     CW_Version.df %>% mutate(Version = 'New Approach CW'))
# #
# # # head(Version.df %>% filter(Version == 'Old Approach'))
# # # head(Version.df %>% filter(Version == 'New Approach JS'))
# # # head(Version.df %>% filter(Version == 'New Approach CW'))
# # # # Plots to compare percentage
# # # png(filename = '~/Desktop/Code/Compare.png')
# # # Version.df %>% ggplot(aes(x = Scenario, y = Frequency)) +
# # #   geom_line(aes(colour = Policy)) +
# # #   facet_grid(Version~.) +
# # #   theme_minimal()
# # # dev.off()
# #
# # # write.csv(file = '~/Desktop/Code/Version.csv',Version.df)
# #
# #
# #
# #
# #
# # # ## Figures
# # # This is the most concise figure script, created 12/28/16 by JS
# # # Makes all five primary figures for MSP Aquaculture Paper
# # # Uses maps created by Becca Gentry in GIS
# # current.directory.data <- paste0(model_directory,'Data')
# # ## Figure dimensions
# # width=7
# # height=5.5
# # res=2400
# # units='in'
# # ## Set Directory
# # cols <- c("Mussel"="Blue","Finfish"="Salmon","Kelp"="Green","Halibut"="Burlywood","Viewshed"='Cyan',"Benthic"='Orange',"Disease"='Black')
# # text.size<-12
# # patch.size=1.5
# # aMatrix <- read.csv(file='~/Desktop/Code/MSP Planning Results April 2016/aMatrix.csv',header=F)
# # # Static.plans.data <- read.csv(file='Static_plans.csv',header=F)
# # Static.values.data <- read.csv(file='~/Desktop/Code/MSP Planning Results April 2016/Dynamic_Values_Export.csv',header=F)
# # colnames(Static.values.data)=c('Mussel','Finfish','Kelp','Halibut','Viewshed','Benthic','Disease')
# # # Static.plans.case.study.data <- read.csv(file='~/Desktop/Code/MSP Planning Results April 2016/April Results/Static_plans_case_study.csv',header=F)
# # # Static.values.case.study.data <- read.csv(file='~/Desktop/Code/MSP Planning Results April 2016/Static_values_case_study.csv',header=F)
# # # Static.percentage.case.study.data <- read.csv('~/Desktop/Code/MSP Planning Results April 2016/Static_percentage_case_study.csv',header=F)
# # # colnames(Static.values.data)=c('Mussel','Finfish','Kelp','Halibut','Viewshed','Benthic','Disease')
# # # colnames(Static.values.case.study.data)=c('Mussel','Finfish','Kelp','Halibut','Viewshed','Benthic','Disease')
# # # data.MSP = list(aMatrix = aMatrix,Static.plans = Static.plans.data,Static.plans.case.study = Static.plans.case.study.data,
# # #               Static.values.case.study = Static.values.case.study.data,Static.percentage.case.study = Static.percentage.case.study.data)
# # # Case Study Stuff
# # I<-which(Static.values.data$Mussel>=.05&Static.values.data$Finfish>=.05&Static.values.data$Kelp>=.05&
# #                            Static.values.data$Halibut>=.95&Static.values.data$Viewshed>=.95&Static.values.data$Benthic>=.95&
# #                            Static.values.data$Disease>=.95)
# # df <- read.csv(file='~/Desktop/Code/SeaGrant_data_complete_2015.csv',stringsAsFactors=FALSE)
# # V1 <- read.csv(file='~/Desktop/Code/V1.csv',header=FALSE)
# # names(df)
# # paste('Cheapest Mussel Farm Costs',df %>% filter(V1 == 1) %>%
# #   select(mussel.annual.operating.costs) %>%
# #   filter(mussel.annual.operating.costs > 0) %>% arrange(mussel.annual.operating.costs) %>% filter(mussel.annual.operating.costs == min(mussel.annual.operating.costs)),'annually')
# #
# # # Case.Study.August <- read.csv(file='Case Study August.csv',header=F)
# # # Static.values.data[I[576],]
# # # summary(factor(Static.plans.data[,I[576]]))
# # # Conventional Models
# # Unconstrained_data <- read.csv(file='~/Desktop/Code/MSP Planning Results April 2016/Unconstrained_Dynamic_Values_April.csv',header=F)
# # names(Unconstrained_data)=c('Mussel','Finfish','Kelp','Halibut','Viewshed','Benthic','Disease')
# # Constrained_data <- read.csv(file='~/Desktop/Code/MSP Planning Results April 2016/Constrained_Dynamic_Values_April.csv',header=F)
# # names(Constrained_data)=c('Mussel','Finfish','Kelp','Halibut','Viewshed','Benthic','Disease')
# #
# # # Unconstrained_data <- read.csv(file='Unconstrained_Static_Values_April.csv',header=F)
# # # names(Unconstrained_data)=c('Mussel','Finfish','Kelp','Halibut','Viewshed','Benthic','Disease')
# # # Constrained_data <- read.csv(file='Constrained_Static_Values_April.csv',header=F)
# # # names(Constrained_data)=c('Mussel','Finfish','Kelp','Halibut','Viewshed','Benthic','Disease')
# #
# # # Pure Profit Conventional Models
# # # Unconstrained.Pure.Profit.data <- read.csv(file='Pure_Profit_Unconstrained_Dynamic_Values_April.csv',header=F)
# # # names(Unconstrained.Pure.Profit.data)=c('Mussel','Finfish','Kelp','Halibut','Viewshed','Benthic','Disease')
# # # Constrained.Pure.Profit.data <- read.csv(file='Pure_Profit_Constrained_Dynamic_Values_April.csv',header=F)
# # # names(Constrained.Pure.Profit.data)=c('Mussel','Finfish','Kelp','Halibut','Viewshed','Benthic','Disease')
# # Master.matrix.max=rbind(Static.values.data,Unconstrained_data,Constrained_data)
# # names(Master.matrix.max)<-c('Mussel','Finfish','Kelp','Halibut','Viewshed','Benthic','Disease')
# # color.vector.max=c(rep('coral',length.out=nrow(Static.values.data)),rep('purple',ts=nrow(Unconstrained_data)),rep('green',ts=nrow(Constrained_data)))
# #
# # # Calculate new annuity values
# # raw_values.df <- read.csv('~/Desktop/CrowTOv1/Raw_Patch_Data.csv')
# # names(raw_values.df) <- c('M','F','K','H','V_F','V_MK','B','D')
# # r <- 0.05
# # t <- 10
# # sector_totals.df <- apply(raw_values.df %>% select(M,F,K,H) %>%
# #   mutate(M = (r*M)/(1-((1+r)^-t)),
# #                 F = (r*F)/(1-((1+r)^-t)),
# #                   K = (r*K)/(1-((1+r)^-t)),
# #                   H = (r*H)/(1-((1+r)^-t))),MARGIN=2,FUN = sum)
# # raw.annuities <- raw_values.df %>% select(M,F,K,H,V_F,V_MK,B,D) %>%
# #   mutate(M = (r*M)/(1-((1+r)^-t)),
# #                 F = (r*F)/(1-((1+r)^-t)),
# #                   K = (r*K)/(1-((1+r)^-t)),
# #                   H = (r*H)/(1-((1+r)^-t)))
# # M_group <- gsub(paste(c('\\(','\\]'),collapse = '|'),'',
# #   names(split(raw.annuities$M,cut(raw.annuities$M,seq(min(raw.annuities$M[which(raw.annuities$M > 0)]),
# #     max(raw.annuities$M),length.out=9)))))
# # F_group <- gsub(paste(c('\\(','\\]'),collapse = '|'),'',
# #   names(split(raw.annuities$F,cut(raw.annuities$F,seq(min(raw.annuities$F[which(raw.annuities$F > 0)]),
# #     max(raw.annuities$F),length.out=9)))))
# # K_group <- gsub(paste(c('\\(','\\]'),collapse = '|'),'',
# #   names(split(raw.annuities$K,cut(raw.annuities$K,seq(min(raw.annuities$K[which(raw.annuities$K > 0)]),
# #     max(raw.annuities$K),length.out=9)))))
# # #### Plot Data
# # current.directory.scripts='~/Desktop/Code/MS Figures/'
# # setwd(current.directory.scripts)
# # theme = theme(plot.margin = unit(c(.2,.2,.2,.2), units = "lines"),
# #                 axis.text = element_blank(),
# #                 axis.title = element_blank(),
# #                 axis.ticks = element_blank(),
# #                 axis.ticks.length = unit(0, "lines"),
# #                 axis.ticks.margin = unit(0, "lines"),
# #                 panel.background=element_rect(fill="white"),
# #                 panel.grid=element_blank(),
# #                 plot.title=element_text(hjust=0))
# #     labs = labs(x = NULL, y = NULL)
# # # png(paste0(current.directory.scripts,'MS_Figures.png'),onefile = T,units=units,width=width, height=height, res=res)
# # for(itor in 1:2){
# # if(itor==1){
# #   png(paste0(current.directory.scripts,'PNGs/Fig 1.png'),units=units,width=width, height=height, res=res)
# # }else{
# #   pdf(paste0(current.directory.scripts,'PDFs/Fig 1.pdf'), width=width, height=height,paper='legal')
# # }
# #   # img <- readTIFF("fig1_Stevens_v3.tif",native=T,info=T)
# #   # g <- rasterGrob(img, interpolate=TRUE)
# #   img <- readPNG("~/Desktop/Code/Fig1A Capture.png",native=T,info=T)
# #   g <- rasterGrob(img, interpolate=TRUE)
# #
# #   a<-qplot(1:10, 1:10, geom="blank") + annotation_custom(g, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) + ggtitle('A') +
# #     theme + labs
# #
# #   img_mussel <- readPNG("~/Desktop/Code/MusselValueApril.png",native=T,info=T)
# #   g_mussel <- rasterGrob(img_mussel, interpolate=TRUE)
# #
# #   foo <- 5.58
# #   store <- NULL
# #   for(itor in 1:9){
# #     foo <- foo - .341
# #     store[itor] <- foo
# #     # print(a)
# #   }
# #
# #   b <- qplot(1:10, 1:10, geom="blank") + ggtitle('B') + annotation_custom(g_mussel, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +
# #     theme +
# #     # annotate(geom='text',x = c(1.75,1.83,1.83,1.83,1.92,1.92,1.92,1.92,1.89),
# #     #   y = store,label=c('0-6','6-6.5','6.5-7','7-7.5','7.5-8.0','8.0-8.5','8.5-9.0','9.0-9.5','9.5-10')) +
# #     labs
# #   # ggsave(filename = "~/Desktop/Code/MusselValueApril_Edit.png",b,device = 'png',dpi = res)
# #
# #   foo <- 5.57
# #   store <- NULL
# #   for(itor in 1:9){
# #     foo <- foo - .35
# #     store[itor] <- foo
# #     # print(a)
# #   }
# #
# #   img_finfish <- readPNG("~/Desktop/Code/FishValueApril.png",native=T,info=T)
# #   g_finfish <- rasterGrob(img_finfish, interpolate=TRUE,just='center')
# #
# #   c<-qplot(1:10, 1:10, geom="blank") + ggtitle('C') + annotation_custom(g_finfish, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +
# #     theme +
# #     # annotate(geom='text',x = c(1.87,1.97,1.92,1.92,1.92,1.92,1.92,1.92,1.92),
# #     #   y = store,label=c('0-0.01','0.01-0.4','0.4-0.8','0.8-1.2','1.2-1.6','1.6-2.0','2.0-2.4','2.4-2.8','2.8-3.2')) +
# #     labs
# #   # ggsave(filename = "~/Desktop/Code/FishValueApril_Edit.png",c,device = 'png',dpi = res)
# #
# #   img_kelp <- readPNG("~/Desktop/Code/KelpValueApril.png",native=T,info=T)
# #   g_kelp <- rasterGrob(img_kelp, interpolate=TRUE, just='center')
# #
# #   # foo <- 5
# #   # store <- NULL
# #   # for(itor in 1:9){
# #   #   foo <- foo - 0.375
# #   #   store[itor] <- foo
# #   #   # print(a)
# #   # }
# #
# #   foo <- 5.58
# #   store <- NULL
# #   for(itor in 1:10){
# #     foo <- foo - .315
# #     store[itor] <- foo
# #     # print(a)
# #   }
# # # annotate(geom='text',x = c(rep(.75,ts=7),.75+.09,.75+.09)
# # # store[length(store) - 2] <- store[length(store) - 2] + .1
# # # store[length(store) - 2] <- store[length(store) - 1] + .2
# #   d<-qplot(1:10, 1:10, geom="blank") + ggtitle('D') + annotation_custom(g_kelp, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +
# #     theme +
# #     # annotate(geom='text',x = c(1.75,1.83,1.83,1.83,1.92,1.92,1.92,1.92,1.89),
# #     #   y = store,label=c('0-3','3-4','4-5','5-6','6-7','7-8','8-9','9-10','10-12'),size = 3.85) +
# #      # annotate(geom='text',x = c(1.75,1.75,1.75,1.75,1.75,1.75,1.75,1.75,1.75,1.80),
# #      #  y = store,label=c('0-1','1-2','2-3','3-4','4-5','5-6','6-7','7-8','8-9','10-12')) +
# #     labs #+
# #     # scale_x_continuous(breaks = seq(0, 10, .25)) +
# #     # scale_y_continuous(breaks = seq(0, 10, .25)) +
# #     # theme(panel.ontop=TRUE,panel.background = element_rect(colour = NA,fill="transparent"),panel.grid.minor = element_blank(),panel.grid.major=element_line(color = 'black'))
# #   # ggsave(filename = "~/Desktop/Code/KelpValueApril_Edit.png",d)
# #
# #   grid.arrange(a,b,c,d,ncol=2,nrow=2)
# #                  # bottom = textGrob(expression(bold('Fig 1: ')~plain('Study domain, spatial constraints and potential value for aquaculture development. (A) Spatial constraints to aquaculture development in the Southern California Bight. (B-D) Potential value (10- year NPV) in each developable cell for mussel, finfish, and kelp aquaculture sectors.')),
# #                  #                   x=1,just='left'))
# #   if(itor == 1){
# #       dev.off()
# #     }else{
# #       dev.off()
# #     }
# # }
# # # Figure 2
# # # pdf(paste0(current.directory.scripts,'Fig 2.pdf'),width=8, height=6.4,paper='legal')
# # # img_MSP <- readPNG(paste0(current.directory.scripts,'Fig 2.png'),native=T,info=T)
# # #   g_MSP <- rasterGrob(img_MSP, interpolate=TRUE, just='center')
# # #   qplot(1:10, 1:10, geom="blank") + annotation_custom(g_MSP, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +
# # #     theme + labs
# # # dev.off()
# #   # png(paste0(current.directory.scripts,'MS_Figures.png'),onefile = T, width=8, height=6.4,res=res,units=units)
# # for(itor in 1){
# #   if(itor==1){
# #     png(paste0(current.directory.scripts,'PNGs/Fig 2.png'),width=8, height=6.4,res=res,units=units)
# #   }else{
# #     pdf(paste0(current.directory.scripts,'PDFs/Fig 2.pdf'),width=8, height=6.4,paper='legal')
# #   }
# #   panel.EF<-function (x, y, itor=0, epsilon=.001, bg = NA, pch = 20, cex = .01, ...)
# #   {
# #     x.MSP=x[color.vector=='coral']
# #     y.MSP=y[color.vector=='coral']
# #
# #     points(x.MSP,y.MSP, pch = 16, col = alpha("dodgerblue",1/75),cex = cex/2)
# #     x.U=x[color.vector=='purple']
# #     y.U=y[color.vector=='purple']
# #     x.S=x[color.vector=='green']
# #     y.S=y[color.vector=='green']
# #     x.EF=NULL
# #     y.EF=NULL
# #     alpha.mat.tmp=seq(from=0,by=epsilon,to=1)
# #     # MSP
# #     for(itor in 1:length(alpha.mat.tmp)){
# #       alpha.tmp=alpha.mat.tmp[itor]
# #       A=(alpha.tmp*x.MSP)+((1-alpha.tmp)*y.MSP)
# #       I=which(A==max(A))
# #       x.EF[itor]=max(unique(x.MSP[I]))
# #       I.tmp.x=which(x.MSP==max(unique(x.MSP[I])))
# #       I.tmp.y=which(y.MSP[I.tmp.x]==max(unique(y.MSP[I.tmp.x])))
# #       y.EF[itor]=unique(y.MSP[I.tmp.x[I.tmp.y]])}
# #     x.EF.original=x.EF;y.EF.original=y.EF;
# #     if(length(unique(x.EF.original))!=1&length(unique(x.EF.original))!=1){
# #       EF.inter=approx(x.EF.original,y=y.EF.original,n=length(alpha.mat.tmp))
# #       x.EF=EF.inter$x;y.EF=EF.inter$y;
# #     }else{
# #     }
# #     lines(sort(x.EF),y.EF[order(x.EF)],col="midnightblue",lwd=2,lty=1)
# #     lines(sort(x.U),y.U[order(x.U)],col = "mediumorchid1",lwd=2,lty=1)
# #     lines(sort(x.S),y.S[order(x.S)],col = "coral1",lwd=2,lty=1)
# #   }
# #   #   pdf.options(width = 8, height = 6.4)
# #   source('~/Desktop/Code/pairs2.R')
# #   color.vector=color.vector.max
# #    # Color Vector For Seperating the MSP from Conventional Solutions
# #   # sample <- rbind(MM_test.df %>% filter(Set == 'MSP') %>% sample_n(size = 1000),
# #   #   MM_test.df %>% filter(Set == 'U') %>% sample_n(size = 500),
# #   #   MM_test.df %>% filter(Set == 'C') %>% sample_n(size = 500))
# #   pairs2(100*Master.matrix.max,lower.panel=panel.EF,
# #          upper.panel=NULL,col=color.vector,cex=0.8,xlim=c(0,100),
# #          ylim=c(0,100),pch=16,font.labels=3,cex.axis=1,las=1,xaxp=c(0,100,4),yaxp=c(0,100,2),
# #          gap=1)
# #   # title(xlab='% of Maximum',line = 1)
# #   title(ylab='% of Maximum')
>>>>>>> e9a6acd68e4a0fd08b18829858208289c9350661
# #   par(xpd=T)
# #   l1<-legend(.33,1,
# #              legend=c('7D Frontier','2D Frontier'),fill=c("lightblue1","midnightblue"),
# #              cex=.75,title=expression(bold('Marine Spatial Planning (MSP)')),
# #              title.adj = 0, bty = 'n', adj = 0, text.width=.25)
# #   l2<-legend(x = l1$rect$left+.0020, y = with(l1$rect, top - h)-.005,
# #              legend=c('Constrained','Unconstrained'),fill=c("coral1","mediumorchid1"),
# #              cex=.75,title=expression(bold('Conventional Planning ')),
# #              title.adj = 0, bty = 'n', adj = 0, text.width=.25)
<<<<<<< HEAD
# #   title(xlab='[%] of Maximum',line=1)
# #   title(ylab='[%] of Maximum')
# # dev.off()#}
# # beep()
#
# # source('figure_S28_code.R')
=======
# #   inset.figure.proportion = 1/3
# #   inset.figure.dims = c(rep(width*(inset.figure.proportion),ts = 2))
# #   subplot(source(file='~/Desktop/Code/Tradeoff Cartoon.R'),x='topright',size = inset.figure.dims, type='plt', par = list(cex.main=2.5, cex = .45, lwd = 1))
# #   par(oma=c(0,2,2,0))
# #   title('A', adj = 0, outer = T, cex = .75)
# #   title(xlab='% of Maximum',line = 3.5)
# #   if(itor == 1){
# #     dev.off()
# #   }else{
# #     dev.off()
# #   }
# # }
# # # Figure 3
# # for(itor in 1:2){
# #   if(itor==1){
# #     png(paste0(current.directory.scripts,'PNGs/Fig 3.png'),width=width, height=height,res=res,units=units)
# #   }else{
# #     pdf(paste0(current.directory.scripts,'PDFs/Fig 3.pdf'),width=width, height=height,paper='legal')
# #   }
# # # png(paste0(current.directory.scripts,'MS_Figures.png'),onefile = T, width=width, height=height,res=res,units=units)
# #
# #   img_hot_all <- readPNG("~/Desktop/Code/HotSpots_ALL_Joel.png",native=T,info=T)
# #   g_hot_all <- rasterGrob(img_hot_all, interpolate=TRUE)
# #
# #   a<-qplot(1:10, 1:10, geom="blank") + ggtitle('A') +
# #     annotation_custom(g_hot_all, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)+
# #     theme + labs
# #
# #   img_hot_mussel <- readPNG("~/Desktop/Code/HotSpots_Mussels_Joel.png",native=T,info=T)
# #   g_hot_mussel <- rasterGrob(img_hot_mussel, interpolate=TRUE)
# #
# #   b<-qplot(1:10, 1:10, geom="blank") + ggtitle('B') +
# #     annotation_custom(g_hot_mussel, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)+
# #     theme + labs
# #
# #   img_hot_finfish <- readPNG("~/Desktop/Code/HotSpots_Fish_Joel.png",native=T,info=T)
# #   g_hot_finfish <- rasterGrob(img_hot_finfish, interpolate=TRUE,just='center')
# #
# #   c<-qplot(1:10, 1:10, geom="blank") + ggtitle('C') +
# #     annotation_custom(g_hot_finfish , xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)+
# #     theme + labs
# #
# #   img_hot_kelp <- readPNG("~/Desktop/Code/HotSpots_Kelp_Joel.png",native=T,info=T)
# #   g_hot_kelp <- rasterGrob(img_hot_kelp, interpolate=TRUE, just='center')
# #
# #   d<-qplot(1:10, 1:10, geom="blank") + ggtitle('D') +
# #     annotation_custom(g_hot_kelp, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)+
# #     theme + labs
# #   grid.arrange(a,b,c,d,ncol=2,nrow=2,padding=unit(.1,'line'))
# #   dev.off()
# # }
# # #  pdf(paste0(current.directory.scripts,'Fig 4.pdf'),width=width, height=height,paper='legal')
# # # img_MSP <- readPNG(paste0(current.directory.scripts,'PNGs/Fig 4.png'),native=T,info=T)
# # # g_MSP <- rasterGrob(img_MSP, interpolate=TRUE, just='center')
# # # qplot(1:10, 1:10, geom="blank") + annotation_custom(g_MSP, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +
# # #   theme + labs
# # # dev.off()
# #
# # # Figure 4
# # for(itor in 1:2){
# #   if(itor==1){
# #     png(paste0(current.directory.scripts,'PNGs/Fig 4.png'),width=width, height=height,res=res,units=units)
# #   }else{
# #     pdf(paste0(current.directory.scripts,'PDFs/Fig 4.pdf'),width=width, height=height,paper='legal')
# #   }
# #   Low.impact.solutions=Static.values.data[I,]
# #   LI.names=names(Low.impact.solutions)
# #   ID=seq(to=nrow(Low.impact.solutions),from=1,by=1)
# #   size.tmp=dim(Low.impact.solutions)
# #   Mussel.LI=Low.impact.solutions[,1]
# #   Finfish.LI=Low.impact.solutions[,2]
# #   Kelp.LI=Low.impact.solutions[,3]
# #   Halibut.LI=Low.impact.solutions[,4]
# #   View.LI=Low.impact.solutions[,5]
# #   Benthic.LI=Low.impact.solutions[,6]
# #   Disease.LI=Low.impact.solutions[,7]
# #   Sector.LI=c(rep(LI.names[1],ts=size.tmp[1]),
# #               rep(LI.names[2],ts=size.tmp[1]),
# #               rep(LI.names[3],ts=size.tmp[1]),
# #               rep(LI.names[4],ts=size.tmp[1]),
# #               rep(LI.names[5],ts=size.tmp[1]),
# #               rep(LI.names[6],ts=size.tmp[1]),
# #               rep(LI.names[7],ts=size.tmp[1]))
# #   for(itor in 1:ncol(Low.impact.solutions)){
# #     ID.LI.tmp=1:size.tmp[1]
# #     if(itor>1){
# #       ID.LI=c(ID.LI,ID.LI.tmp)
# #     }else{
# #       ID.LI=ID.LI.tmp
# #     }
# #   }
# #   value.LI.tmp=c(Mussel.LI,Finfish.LI,Kelp.LI,Halibut.LI,View.LI,Benthic.LI,Disease.LI)
# #   Low.impact.solutions=data.frame(ID.LI,value.LI.tmp,Sector.LI)
# #   names(Low.impact.solutions)=c('ID','Value','Sector')
# #   Low.impact.solutions$Sector=factor(Low.impact.solutions$Sector, levels=c('Mussel','Finfish','Kelp','Halibut','Viewshed','Benthic','Disease'))
# #   p.bar<-ggplot(data = Low.impact.solutions,aes(x=ID,y=Value,fill=Sector,color=Sector))+
# #     geom_bar(stat="identity")+facet_grid(.~Sector)+ggtitle('A')+
# #     scale_y_continuous(labels = percent,limits=c(0,1))+
# #     scale_fill_manual(values=cols)+
# #     scale_color_manual(values=cols)+
# #     theme(axis.title.x=element_blank(),
# #           axis.text.x=element_blank(),
# #           axis.ticks.x=element_blank(),
# #           panel.background=element_rect(color='white',fill='white'),panel.spacing = unit(.75, "lines"),
# #           strip.background=element_rect(fill='white'),strip.text=element_text(size=text.size*.75),
# #           panel.grid=element_blank(),axis.title.y=element_text(size=text.size,color="black"),
# #           axis.text.y=element_text(size=text.size,color="black"),legend.title=element_text(size=text.size*.75),
# #           legend.text=element_text(size=text.size),plot.title=element_text(hjust=0))
# #
# #   img_case_study_percent <- readPNG("~/Desktop/Code/CaseStudyPercentage_August.png",native=T,info=T)
# #   g_case_study_percent <- rasterGrob(img_case_study_percent, interpolate=TRUE,just='center')
# #
# #   p.map<-qplot(1:10, 1:10, geom="blank") + ggtitle('B') +
# #     annotation_custom(g_case_study_percent, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)+
# #     theme(plot.margin = unit(c(.2,.2,.4,.2), units = "lines"),axis.text.x = element_blank(),axis.text.y = element_blank(),axis.ticks = element_blank(),
# #           axis.title=element_blank(),panel.background=element_rect(fill="white"),panel.grid=element_blank(),
# #           plot.title=element_text(hjust=0))
# #
# #   img_case_study_plan <- readPNG("~/Desktop/Code/CaseStudySpecies_August.png",info=T)
# #   g_case_study_plan <- rasterGrob(img_case_study_plan, interpolate=TRUE,just='center')
# #
# #   p.map.case.study<-qplot(1:10, 1:10, geom="blank") + ggtitle('C') +
# #     annotation_custom(g_case_study_plan , xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)+
# #     theme(plot.margin = unit(c(.2,.2,.4,.2), units = "lines"),axis.text.x = element_blank(),axis.text.y = element_blank(),axis.ticks = element_blank(),
# #           axis.title=element_blank(),panel.background=element_rect(fill="white"),panel.grid=element_blank(),
# #           plot.title=element_text(hjust=0))
# #
# #   p.bar.case.study <- ggplot(data = Low.impact.solutions[Low.impact.solutions$ID==576,],
# #                              aes(x=ID,y=Value,fill=Sector,color=Sector))+ggtitle('D')+
# #     geom_bar(stat="identity")+facet_grid(.~Sector)+
# #     scale_y_continuous(labels = percent,limits=c(0,1))+
# #     scale_fill_manual(values=cols)+
# #     scale_color_manual(values=cols)+
# #     theme(axis.ticks.x=element_blank(),axis.title.x=element_blank(),
# #           axis.title.y=element_blank(),axis.text.x=element_blank(),
# #           panel.background=element_rect(color='white',fill='white'),
# #           strip.background=element_rect(fill='white'),strip.text=element_blank(),
# #           panel.grid=element_blank(),
# #           axis.text.y=element_blank(),axis.ticks=element_blank(),legend.position='none',plot.title=element_text(hjust=0),plot.margin = unit(c(.2,.2,.2,.2), "cm"))
# #   g=ggplotGrob(p.bar.case.study)
# #   p.map.case.study.full<-p.map.case.study+annotation_custom(grob=g,xmin=5.5,xmax=10.25,ymin=7.5)
# #   grid.arrange(p.bar,p.map,p.map.case.study.full,ncol=2,nrow=2,layout_matrix = rbind(c(1,1),c(2,3)))
# #   dev.off()
# # }
# #   # pdf.options(width = 9.5, height = 7)
# # # Figure 5
# # for(itor in 1:2){
# #   if(itor==1){
# #     png(paste0(current.directory.scripts,'PNGs/Fig 5.png'),width=8, height=8,res=res,units=units)
# #   }else{
# #     pdf(paste0(current.directory.scripts,'PDFs/Fig 5.pdf'),width=8, height=8,paper='legal')
# #   }
# #   source('~/Desktop/Code/value.of.MSP.loop_interpol.R')
# #   source('~/Desktop/Code/value.of.MSP.fx_2_interpol.R')
# #   source('~/Desktop/Code/figure_5_code.R')
>>>>>>> e9a6acd68e4a0fd08b18829858208289c9350661
# #   MSP.value.data=rbind(Mussel.value.tmp,Finfish.value.tmp,Kelp.value.tmp)
# #   MSP.value.data.points=rbind(Mussel.value.tmp.points,Finfish.value.tmp.points,Kelp.value.tmp.points)
# #   names(MSP.value.data)=c('Aquaculture.Value','Value.of.MSP','Group','Sector.Name','Sector.Type','Type.of.Conventional')
# #   names(MSP.value.data.points)=c('Aquaculture.Value','Value.of.MSP','Group','Sector.Name','Sector.Type','Type.of.Conventional')
# #   MSP.value.data$Sector.Name=factor(MSP.value.data$Sector.Name,levels=c('Mussel','Finfish','Kelp','Halibut','Viewshed','Benthic','Disease'))
# #   Value.of.MSP.grid.plot<-ggplot(data=MSP.value.data)+
# #     geom_point(data=subset(MSP.value.data[as.integer(MSP.value.data$Group)!=as.integer(MSP.value.data$Sector.Name),]),
# #                aes(x=Aquaculture.Value,y=Value.of.MSP,shape=Type.of.Conventional,color=Type.of.Conventional),size=1.25)+
# #     facet_grid(Sector.Name~Group)+scale_y_continuous(labels = percent,limits=c(0,1),breaks=c(0,.5,1))+
# #     scale_x_continuous(labels = percent,breaks=c(0,.5,1))+
# #     geom_text(data=subset(MSP.value.data[as.integer(MSP.value.data$Group)==as.integer(MSP.value.data$Sector.Name),]),x=.5,y=.5,size=15,label='NA')+
# #     geom_point(data=MSP.value.data.points,aes(x=Aquaculture.Value,y=Value.of.MSP,shape=Type.of.Conventional,color=Type.of.Conventional),size=1.2)+
<<<<<<< HEAD
# #     xlab("Aquaculture Value")+ylab("Value of MSP")+
=======
# #     xlab("Aquaculture Value")+ylab("Value of Marine Spatial Planning")+
>>>>>>> e9a6acd68e4a0fd08b18829858208289c9350661
# #     scale_colour_manual(name = "Conventional Planning :",labels=c("Constrained","Unconstrained"),values=c("coral1","mediumorchid1"))+
# #     scale_shape_manual(name = "Conventional Planning :",labels=c("Constrained","Unconstrained"),values=c(16,24))+
# #     scale_fill_manual(name = "Conventional Planning :",labels=c("Constrained","Unconstrained"),values=c("coral1","mediumorchid1"))+
# #     #     geom_line(aes(x=c(0,1),y=c(0,0)),color="grey",linetype='dashed',size=1)+
# #     theme_bw(base_size = 15)+theme(panel.grid=element_blank(),legend.position="bottom")
<<<<<<< HEAD
# #   Fig_S10<-Value.of.MSP.grid.plot+theme(axis.title=element_text(size=12),panel.spacing = unit(1, "lines"),
# #                                strip.background=element_rect(fill='white'),strip.text=element_text(size=10),axis.ticks.margin=unit(1,'lines'))
#
# #   for(itor in 1:2){
# #     if(itor==1){
# #       png('Fig S10.png',width=8, height=8,res=res,units=units)
# #     }else{
# #       pdf('Fig S10.pdf',width=8, height=8,paper='legal')
# #     }
# #     print(Fig_S10)
# #   dev.off()}
=======
# #   Value.of.MSP.grid.plot+theme(axis.title=element_text(size=12),panel.spacing = unit(1, "lines"),
# #                                strip.background=element_rect(fill='white'),strip.text=element_text(size=10),axis.ticks.margin=unit(1,'lines'))
# #   dev.off()
# # }
# #
# # # ## Combine Figures
# # # require(png)
# # # require(grid)
# # # png(paste0(current.directory.scripts,'/PNGs/MS_Figures.png'))
# # # lapply(ll <- list.files(path = '~/Desktop/Code/MS Figures/PNGs',patt='.*[.]png'),function(x){
# # #   img <- as.raster(readPNG(paste0('~/Desktop/Code/MS Figures/PNGs/',x)))
# # #   grid.newpage()
# # #   grid.raster(img, interpolate = FALSE)
# #
# # # })
# # # dev.off()
# #
# #
# #
# #
# # # plots <- lapply(ll <- list.files(path = paste0(current.directory.scripts,'PNGs/'), patt='.*[.]png'),function(x){
# # #   img <- as.raster(readPNG(x))
# # #   rasterGrob(img, interpolate = FALSE)
# # # })
# # # require(ggplot2)
# # # require(gridExtra)
# # # ggsave("multipage.pdf", marrangeGrob(grobs=plots, nrow=2, ncol=2))
# #
# #
# #
# # ## Figure 1
# # ## SI Figures
# # current.directory.code <- '~/Desktop/Code/SI Figures Code/'
# # current.directory.figures <- '~/Desktop/Code/SI Figures/'
# # # Coastline.data<-read.csv(file='Coastline.csv',header=F)
# # theme_3 = theme(plot.margin = unit(c(.1,.1,.1,.1), units = "lines"),
# #               axis.text = element_blank(),
# #               axis.title = element_blank(),
# #               axis.ticks = element_blank(),
# #               # axis.ticks.length = unit(0, "lines"),
# #               # axis.ticks.margin = unit(0, "lines"),
# #               panel.background=element_rect(fill="white"),
# #               panel.grid=element_blank(),
# #               plot.title=element_text(hjust =.35,vjust=.1))
# #
# # theme_2 = theme(plot.margin = unit(c(.1,.1,.1,.1), units = "lines"),
# #               axis.text = element_blank(),
# #               axis.title = element_blank(),
# #               axis.ticks = element_blank(),
# #               # axis.ticks.length = unit(0, "lines"),
# #               # axis.ticks.margin = unit(0, "lines"),
# #               panel.background=element_rect(fill="white"),
# #               panel.grid=element_blank(),
# #               plot.title=element_text(hjust =.30,vjust=.1))
# #
# # labs = labs(x = NULL, y = NULL)
# # # Figure S1
# # for(itor in 1:2){
# #   if(itor==1){
# #     png(paste0(current.directory.figures,'PNGs/Fig S1.png'),width=width, height=height,res=res,units=units)
# #   }else{
# #     pdf(paste0(current.directory.figures,'PDFs/Fig S1.pdf'),width=width, height=height,paper='legal')
# #   }
# #   img_s1 <- readPNG(paste0(current.directory.code,'Fig S1.png'),native=T,info=T)
# #   g_s1 <- rasterGrob(img_s1, interpolate=TRUE)
# #
# #   s1<-qplot(1:10, 1:10, geom="blank") +
# #     annotation_custom(g_s1, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)+
# #     theme_3 + labs
# #   print(s1)
# #   dev.off()}
# # # Figure S2
# # for(itor in 1:2){
# #   if(itor==1){
# #     png(paste0(current.directory.figures,'PNGs/Fig S2.png'),width=width, height=height,res=res,units=units)
# #   }else{
# #     pdf(paste0(current.directory.figures,'PDFs/Fig S2.pdf'),width=width, height=height,paper='legal')
# #   }
# #   img_S2A <- readPNG(paste0(current.directory.code,"S2A.png"),native=T,info=T)
# #   g_S2A <- rasterGrob(img_S2A, interpolate=TRUE)
# #
# #   S2A<-qplot(1:10, 1:10, geom="blank") + ggtitle('A') +
# #     annotation_custom(g_S2A, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)+
# #     theme_3 + labs
# #
# #   img_S2B <- readPNG(paste0(current.directory.code,"S2B.png"),native=T,info=T)
# #   g_S2B <- rasterGrob(img_S2B, interpolate=TRUE)
# #
# #   S2B<-qplot(1:10, 1:10, geom="blank") + ggtitle('B') +
# #     annotation_custom(g_S2B, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)+
# #     theme_3 + labs
# #
# #   img_S2C <- readPNG(paste0(current.directory.code,"S2C.png"),native=T,info=T)
# #   g_S2C <- rasterGrob(img_S2C, interpolate=TRUE)
# #
# #   S2C<-qplot(1:10, 1:10, geom="blank") + ggtitle('C') +
# #     annotation_custom(g_S2C, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)+
# #     theme_3 + labs
# #
# #   grid.arrange(S2A,S2B,S2C,ncol=1,nrow=3,padding=unit(-.1,'line'))
# #   dev.off()}
# # # S3
# #   for(itor in 1:2){
# #     if(itor==1){
# #       png(paste0(current.directory.figures,'PNGs/Fig S3.png'),width=width, height=height,res=res,units=units)
# #     }else{
# #       pdf(paste0(current.directory.figures,'PDFs/Fig S3.pdf'),width=width, height=height,paper='legal')
# #     }
# #   img_S3A <- readPNG(paste0(current.directory.code,"S3A.png"),native=T,info=T)
# #   g_S3A <- rasterGrob(img_S3A, interpolate=TRUE)
# #
# #   S3A<-qplot(1:10, 1:10, geom="blank") + ggtitle('A') +
# #     annotation_custom(g_S3A, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)+
# #     theme_3 + labs
# #
# #   img_S3B <- readPNG(paste0(current.directory.code,"S3B.png"),native=T,info=T)
# #   g_S3B <- rasterGrob(img_S3B, interpolate=TRUE)
# #
# #   S3B<-qplot(1:10, 1:10, geom="blank") + ggtitle('B') +
# #     annotation_custom(g_S3B, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)+
# #     theme_3 + labs
# #
# #   img_S3C <- readPNG(paste0(current.directory.code,"S3C.png"),native=T,info=T)
# #   g_S3C <- rasterGrob(img_S3C, interpolate=TRUE)
# #
# #   S3C<-qplot(1:10, 1:10, geom="blank") + ggtitle('C') +
# #     annotation_custom(g_S3C, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)+
# #     theme_3 + labs
# #   grid.arrange(S3A,S3B,S3C,ncol=1,nrow=3,padding=unit(-.1,'line'))
# #   dev.off()}
# # # S4
# # for(itor in 1:2){
# #   if(itor==1){
# #     png(paste0(current.directory.figures,'PNGs/Fig S4.png'),width=width, height=height,res=res,units=units)
# #   }else{
# #     pdf(paste0(current.directory.figures,'PDFs/Fig S4.pdf'),width=width, height=height,paper='legal')
# #   }
# #   img_s4 <- readPNG(paste0(current.directory.code,"Fig S4.png"),native=T,info=T)
# #   g_s4 <- rasterGrob(img_s4, interpolate=TRUE)
# #
# #   s4<-qplot(1:10, 1:10, geom="blank") +
# #     annotation_custom(g_s4, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)+
# #     theme_3 + labs
# #   print(s4)
# #   dev.off()}
# # # S5
# #   for(itor in 1:2){
# #     if(itor==1){
# #       png(paste0(current.directory.figures,'PNGs/Fig S5.png'),width=width, height=height,res=res,units=units)
# #     }else{
# #       pdf(paste0(current.directory.figures,'PDFs/Fig S5.pdf'),width=width, height=height,paper='legal')
# #     }
# #
# #   img_S5A <- readPNG(paste0(current.directory.code,"S5A.png"),native=T,info=T)
# #   g_S5A <- rasterGrob(img_S5A, interpolate=TRUE)
# #
# #   S5A<-qplot(1:10, 1:10, geom="blank") + ggtitle('A') +
# #     annotation_custom(g_S5A, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)+
# #     theme_2 + labs
# #
# #   img_S5B <- readPNG(paste0(current.directory.code,"S5B.png"),native=T,info=T)
# #   g_S5B <- rasterGrob(img_S5B, interpolate=TRUE)
# #
# #   S5B<-qplot(1:10, 1:10, geom="blank") + ggtitle('B') +
# #     annotation_custom(g_S5B, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)+
# #     theme_2 + labs
# #
# #   grid.arrange(S5A,S5B,ncol=1,nrow=2,padding=unit(-.1,'line'))
# #   dev.off()}
# # # S6
# # for(itor in 1:2){
# #   if(itor==1){
# #      png(paste0(current.directory.figures,'PNGs/Fig S6.png'),width=width, height=height,res=res,units=units)
# #   }else{
# #      pdf(paste0(current.directory.figures,'PDFs/Fig S6.pdf'),width=width, height=height,paper='legal')
# #   }
# #   img_S6A <- readPNG(paste0(current.directory.code,"S6A.png"),native=T,info=T)
# #   g_S6A <- rasterGrob(img_S6A, interpolate=TRUE)
# #
# #   S6A<-qplot(1:10, 1:10, geom="blank") + ggtitle('A') +
# #     annotation_custom(g_S6A, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)+
# #     theme_2 + labs
# #
# #   img_S6B <- readPNG(paste0(current.directory.code,"S6B.png"),native=T,info=T)
# #   g_S6B <- rasterGrob(img_S6B, interpolate=TRUE)
# #
# #   S6B<-qplot(1:10, 1:10, geom="blank") + ggtitle('B') +
# #     annotation_custom(g_S6B, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)+
# #     theme_2 + labs
# #
# #   grid.arrange(S6A,S6B,ncol=1,nrow=2,padding=unit(-.1,'line'))
# #   dev.off()
# # }
# # # S7
# # for(itor in 1:2){
# #   if(itor==1){
# #     png(paste0(current.directory.figures,'PNGs/Fig S7.png'),width=width, height=height,res=res,units=units)
# #   }else{
# #     pdf(paste0(current.directory.figures,'PDFs/Fig S7.pdf'),width=width, height=height,paper='legal')
# #   }
# #   img_S7A <- readPNG(paste0(current.directory.code,"S7A.png"),native=T,info=T)
# #   g_S7A <- rasterGrob(img_S7A, interpolate=TRUE)
# #
# #   S7A<-qplot(1:10, 1:10, geom="blank") + ggtitle('A') +
# #     annotation_custom(g_S7A, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)+
# #     theme_3 + labs + theme(plot.title=element_text(hjust=0))
# #
# #   img_S7B <- readPNG(paste0(current.directory.code,"S7B.png"),native=T,info=T)
# #   g_S7B <- rasterGrob(img_S7B, interpolate=TRUE)
# #
# #   S7B<-qplot(1:10, 1:10, geom="blank") + ggtitle('B') +
# #     annotation_custom(g_S7B, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)+
# #     theme_3 + labs + theme(plot.title=element_text(hjust=0))
# #
# #   img_S7C <- readPNG(paste0(current.directory.code,"S7C.png"),native=T,info=T)
# #   g_S7C <- rasterGrob(img_S7C, interpolate=TRUE)
# #
# #   S7C<-qplot(1:10, 1:10, geom="blank") + ggtitle('C') +
# #     annotation_custom(g_S7C, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)+
# #     theme_3 + labs + theme(plot.title=element_text(hjust=0))
# #
# #   img_S7D <- readPNG(paste0(current.directory.code,"S7D.png"),native=T,info=T)
# #   g_S7D <- rasterGrob(img_S7D, interpolate=TRUE)
# #
# #   S7D<-qplot(1:10, 1:10, geom="blank") + ggtitle('D') +
# #     annotation_custom(g_S7D, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)+
# #     theme_3 + labs + theme(plot.title=element_text(hjust=0))
# #
# #   grid.arrange(S7A,S7B,S7C,S7D,ncol=2,nrow=2,padding=unit(.5,'line'))
# #   dev.off()
# # }
# # # S8
# # for(itor in 1:2){
# #   if(itor==1){
# #     png(paste0(current.directory.figures,'PNGs/Fig S8.png'),width=width, height=height,res=res,units=units)
# #   }else{
# #     pdf(paste0(current.directory.figures,'PDFs/Fig S8.pdf'),width=width, height=height,paper='legal')
# #   }
# #   img_S8A <- readPNG(paste0(current.directory.code,"S8A.png"),native=T,info=T)
# #   g_S8A <- rasterGrob(img_S8A, interpolate=TRUE)
# #
# #   S8A<-qplot(1:10, 1:10, geom="blank") + ggtitle('A') +
# #     annotation_custom(g_S8A, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)+
# #     theme_3 + labs
# #
# #   img_S8B <- readPNG(paste0(current.directory.code,"S8B.png"),native=T,info=T)
# #   g_S8B <- rasterGrob(img_S8B, interpolate=TRUE)
# #
# #   S8B<-qplot(1:10, 1:10, geom="blank") + ggtitle('B') +
# #     annotation_custom(g_S8B, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)+
# #     theme_3 + labs
# #
# #   grid.arrange(S8A,S8B,ncol=1,nrow=2,padding=unit(.5,'line'))
# #   dev.off()
# # }
# #
# # # # setwd('C:/Users/Joel/Desktop/Thesis YTD/Code')
# # #   # panel.EF<-function (x, y, itor=0, epsilon=.001, bg = NA, pch = 20, cex = .01, ...)
# # #   # {
# # #   #   x.MSP=x[color.vector=='coral']
# # #   #   y.MSP=y[color.vector=='coral']
# # #   #
# # #   #   points(x.MSP,y.MSP, pch = 16, col = alpha("dodgerblue",1/75),cex = cex/2)
# # #   #   x.U=x[color.vector=='purple']
# # #   #   y.U=y[color.vector=='purple']
# # #   #   x.S=x[color.vector=='green']
# # #   #   y.S=y[color.vector=='green']
# # #   #   x.EF=NULL
# # #   #   y.EF=NULL
# # #   #   alpha.mat.tmp=seq(from=0,by=epsilon,to=1)
# # #   #   # MSP
# # #   #   for(itor in 1:length(alpha.mat.tmp)){
# # #   #     alpha.tmp=alpha.mat.tmp[itor]
# # #   #     A=(alpha.tmp*x.MSP)+((1-alpha.tmp)*y.MSP)
# # #   #     I=which(A==max(A))
# # #   #     x.EF[itor]=max(unique(x.MSP[I]))
# # #   #     I.tmp.x=which(x.MSP==max(unique(x.MSP[I])))
# # #   #     I.tmp.y=which(y.MSP[I.tmp.x]==max(unique(y.MSP[I.tmp.x])))
# # #   #     y.EF[itor]=unique(y.MSP[I.tmp.x[I.tmp.y]])}
# # #   #   x.EF.original=x.EF;y.EF.original=y.EF;
# # #   #   if(length(unique(x.EF.original))!=1&length(unique(x.EF.original))!=1){
# # #   #     EF.inter=approx(x.EF.original,y=y.EF.original,n=length(alpha.mat.tmp))
# # #   #     x.EF=EF.inter$x;y.EF=EF.inter$y;
# # #   #   }else{
# # #   #   }
# # #   #   lines(sort(x.EF),y.EF[order(x.EF)],col="midnightblue",lwd=2,lty=1)
# # #   #   lines(sort(x.U),y.U[order(x.U)],col = "mediumorchid1",lwd=2,lty=1)
# # #   #   lines(sort(x.S),y.S[order(x.S)],col = "coral1",lwd=2,lty=1)
# # #   # }
# # #   #   pdf.options(width = 8, height = 6.4)
# #
# # #   # source('pairs2.R')
# # # # for(itor in 1:2){
# # # #   if(itor==1){
# # #     png('Fig S9.png',width=8, height=6.4,res=res,units=units)
# # #   # }else{
# # #   #   pdf('Fig S9.pdf',width=8, height=6.4,paper='legal')
# # #   # # }
# # #   color.vector=color.vector.max # Color Vector For Seperating the MSP from Conventional Solutions
# # #   pairs2(100*Master.matrix.max.pure.profit,lower.panel=panel.EF,
# # #          upper.panel=NULL,col=color.vector,cex=0.8,xlim=c(0,100),
# # #          ylim=c(0,100),pch=16,font.labels=3,cex.axis=1,las=1,xaxp=c(0,100,4),yaxp=c(0,100,2),
# # #          gap=1)
# # #   par(xpd=T)
# # #   l1<-legend(.33,1,
# # #              legend=c('7D Frontier','2D Frontier'),fill=c("lightblue1","midnightblue"),
# # #              cex=.75,title=expression(bold('Marine Spatial Planning (MSP)')),
# # #              title.adj = 0, bty = 'n', adj = 0, text.width=.25)
# # #   l2<-legend(x = l1$rect$left+.0020, y = with(l1$rect, top - h)-.005,
# # #              legend=c('Constrained','Unconstrained'),fill=c("coral1","mediumorchid1"),
# # #              cex=.75,title=expression(bold('Conventional Planning ')),
# # #              title.adj = 0, bty = 'n', adj = 0, text.width=.25)
# # #   title(xlab='[%] of Maximum',line=1)
# # #   title(ylab='[%] of Maximum')
# # # dev.off()#}
# # # beep()
# #
# # # source('figure_S28_code.R')
# # #   MSP.value.data=rbind(Mussel.value.tmp,Finfish.value.tmp,Kelp.value.tmp)
# # #   MSP.value.data.points=rbind(Mussel.value.tmp.points,Finfish.value.tmp.points,Kelp.value.tmp.points)
# # #   names(MSP.value.data)=c('Aquaculture.Value','Value.of.MSP','Group','Sector.Name','Sector.Type','Type.of.Conventional')
# # #   names(MSP.value.data.points)=c('Aquaculture.Value','Value.of.MSP','Group','Sector.Name','Sector.Type','Type.of.Conventional')
# # #   MSP.value.data$Sector.Name=factor(MSP.value.data$Sector.Name,levels=c('Mussel','Finfish','Kelp','Halibut','Viewshed','Benthic','Disease'))
# # #   Value.of.MSP.grid.plot<-ggplot(data=MSP.value.data)+
# # #     geom_point(data=subset(MSP.value.data[as.integer(MSP.value.data$Group)!=as.integer(MSP.value.data$Sector.Name),]),
# # #                aes(x=Aquaculture.Value,y=Value.of.MSP,shape=Type.of.Conventional,color=Type.of.Conventional),size=1.25)+
# # #     facet_grid(Sector.Name~Group)+scale_y_continuous(labels = percent,limits=c(0,1),breaks=c(0,.5,1))+
# # #     scale_x_continuous(labels = percent,breaks=c(0,.5,1))+
# # #     geom_text(data=subset(MSP.value.data[as.integer(MSP.value.data$Group)==as.integer(MSP.value.data$Sector.Name),]),x=.5,y=.5,size=15,label='NA')+
# # #     geom_point(data=MSP.value.data.points,aes(x=Aquaculture.Value,y=Value.of.MSP,shape=Type.of.Conventional,color=Type.of.Conventional),size=1.2)+
# # #     xlab("Aquaculture Value")+ylab("Value of MSP")+
# # #     scale_colour_manual(name = "Conventional Planning :",labels=c("Constrained","Unconstrained"),values=c("coral1","mediumorchid1"))+
# # #     scale_shape_manual(name = "Conventional Planning :",labels=c("Constrained","Unconstrained"),values=c(16,24))+
# # #     scale_fill_manual(name = "Conventional Planning :",labels=c("Constrained","Unconstrained"),values=c("coral1","mediumorchid1"))+
# # #     #     geom_line(aes(x=c(0,1),y=c(0,0)),color="grey",linetype='dashed',size=1)+
# # #     theme_bw(base_size = 15)+theme(panel.grid=element_blank(),legend.position="bottom")
# # #   Fig_S10<-Value.of.MSP.grid.plot+theme(axis.title=element_text(size=12),panel.spacing = unit(1, "lines"),
# # #                                strip.background=element_rect(fill='white'),strip.text=element_text(size=10),axis.ticks.margin=unit(1,'lines'))
# #
# # #   for(itor in 1:2){
# # #     if(itor==1){
# # #       png('Fig S10.png',width=8, height=8,res=res,units=units)
# # #     }else{
# # #       pdf('Fig S10.pdf',width=8, height=8,paper='legal')
# # #     }
# # #     print(Fig_S10)
# # #   dev.off()}
>>>>>>> e9a6acd68e4a0fd08b18829858208289c9350661
