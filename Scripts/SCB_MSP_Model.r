# //TODO: Integrate call to matlab for running the dynamic halibut results
#//TODO: Make sure the figure scratch file is updated
#//TODO: Get all of the directories finalized
#//TODO: Integrate seeding scripts
# MATLAB requires 'Mapping Toolbox', 'Bioinformatics', 'Parallel Optimization'
# Set current working directory as a string
# Load necessary R libraries. For the function R_Libraries, enter T if this is the first time running the model. This will
# install all of the necessary libraries and load them into the
# current workspace.
# source(paste0(wkdir,'/MSP_Model/Scripts','/R_Libraries.r'))
# R_Libraries(FALSE) # After the first initial run this can be set to F
# Install R markdown
# install.packages("knitr",repos = 'https://cran.mtu.edu/')
# library(knitr)
# Set global variables
n.sector <- 7 # Number of sectors
epsilon <- 0.2 # Stepsize of sector weights
t <- 10 # Time Horizon
r <- 0.05 # Discount rate

# Read sector data
# loadWorkbook(paste0(fdirs$inpdatadir,'SeaGrant_data_complete_2015.csv'))
sector_data.df <- read.csv(paste0(inpdatadir,'SeaGrant_data_complete_2015.csv'),stringsAsFactors = FALSE)

fulldomain <- sector_data.df$TARGET_FID # Model domain
discount_factor <- 1/((1+r)^c(1:t))
# Make a numeric matrix of the discount_factor with the dimensions of
# rows = 6425 (number of sites in the domain), columns = 10 (time horizon)
r_iy_aqua <- do.call(rbind, replicate(length(fulldomain), discount_factor, simplify=FALSE))

# Calculate the 10-year Net Present Value (NPV) and annuities for each form of aquaculture
# and halibut (the only impacted sector with direct monetary value)
  # Function to calculate the NPV/annuities for Mussel and Kelp
  Value.MK <- function(yield,upfront.cost,annual.cost,price,t,r,r_iy_aqua){
    revenue <- do.call(cbind, replicate(t, yield * price, simplify=FALSE))
    cost <- cbind(upfront.cost + annual.cost,
      do.call(cbind, replicate(t - 1, annual.cost, simplify=FALSE)))
    profit <- (revenue - cost) * r_iy_aqua
    NPV <- apply(profit, FUN = sum, MARGIN = 1)
    NPV[NPV < 0] <- 0
    print(length(which(NPV > 0)))
    Annuity = (r*NPV)/(1-((1+r)^-t))
    return(list(NPV = NPV,Annuity = Annuity))
  }
  # Function to calculate the NPV/annuities for Finfish
  Value.F <- function(yield,costs,price,t,r,r_iy_aqua){
    revenue <- do.call(cbind, replicate(t, yield * price, simplify=FALSE))
    cost <- sector_data.df$fish.annual.operating.costs
    profit <- (revenue - cost) * r_iy_aqua
    NPV <- apply(profit, FUN = sum, MARGIN = 1)
    NPV[NPV < 0] <- 0
    print(length(which(NPV > 0)))
    Annuity = (r*NPV)/(1-((1+r)^-t))
    return(list(NPV = NPV,Annuity = Annuity))
  }
# Calculate respective NPV/annuities
# Mussel, fixed price of $3.30 per kg
  M <- Value.MK(sector_data.df$mussel.yield,
    sector_data.df$mussel.upfront.cost,
    sector_data.df$mussel.annual.operating.cost,3.3,t,r,r_iy_aqua)
# Finfish, fixed price of $8.00 per kg
  F <- Value.F(sector_data.df$fish.yield,
    sector_data.df$fish.upfront.cost,
    unique(sector_data.df$fish.price[sector_data.df$fish.price>0]),t,r,r_iy_aqua)
# Kelp, fixed price of $3.00 per kg
  K <- Value.MK(sector_data.df$kelp.yield,
    sector_data.df$kelp.upfront.cost,
    sector_data.df$kelp.annual.operating.cost,3,t,r,r_iy_aqua)
# Remove unprofitable sites and generate seperate vectors for
# Those sites in which ventures will be profitable for at least one type
# of aquaculture (var Aqua.Full.Domain),
Aqua.Full.Domain <- data.frame(M$Annuity,F$Annuity,K$Annuity)
Aqua.Full.Domain.Logical <- apply(1 * (Aqua.Full.Domain > 0),FUN=sum,MARGIN=1) > 0
# Sites that will be profitable for mussel (M.V)
M.V <- M$Annuity[Aqua.Full.Domain.Logical]
# Sites that will be profitable for finfish (F.V)
F.V <- F$Annuity[Aqua.Full.Domain.Logical]
# Sites that will be profitable for kelp (K.V)
K.V <- K$Annuity[Aqua.Full.Domain.Logical]

# Run the Halibut fishing model and then load the results
if(readline("Run halibut model or load results Y/N? ") == 'Y'){
  system2(matlab_root,
    args = c('-nodesktop','-noFigureWindows','-nosplash','-r',
    paste0("run\\(\\'",scrpdir,"Halibut/Tuner_free_params_v9.m\\'\\)")))
}
# sum((r*read.csv(paste0(outdatadir,'Target_FID_and_Yi_fulldomain_NPV_at_MSY_noAqua.csv'),header=FALSE)[,2])/(1-((1+r)^-t)))
H.V  <- (r*read.csv(paste0(outdatadir,'Target_FID_and_Yi_fulldomain_NPV_at_MSY_noAqua.csv'),header=FALSE)[,2][Aqua.Full.Domain.Logical])/(1-((1+r)^-t))
# (r*read.csv(paste0('~/MSP_Model/Input/Data/TargFID_and_YifulldomainNPV_at_MSY_noAqua.csv'),header=FALSE)[,2][Aqua.Full.Domain.Logical])/(1-((1+r)^-t))


# Load Viewshed Data
# sector_data.df$park_views_8k[which(sector_data.df$TARGET_FID %in% viewshed_new$FID)] <- as.character(viewshed_new$NEW8KTOT)
# sector_data.df$park_views_3k[which(sector_data.df$TARGET_FID %in% viewshed_new$FID)] <- as.character(viewshed_new$NEW3KTOT)
# V_F.V <- (as.numeric(gsub(",", "", sector_data.df$res_views_8k)) + as.numeric(gsub(",", "",sector_data.df$park_views_8k)))[Aqua.Full.Domain.Logical]
# V_MK.V  <- (as.numeric(gsub(",", "", sector_data.df$res_views_3k)) + as.numeric(gsub(",", "",sector_data.df$park_view_3k)))[Aqua.Full.Domain.Logical]
viewshed_new <- read.csv(paste0(inpdatadir,'viewshed_full_updated_2017_w_new_rec_views.csv'))
V_F.V <- viewshed_new[viewshed_new$FID %in% fulldomain[Aqua.Full.Domain.Logical],]$NEW8KTOT
V_MK.V <- viewshed_new[viewshed_new$FID %in% fulldomain[Aqua.Full.Domain.Logical],]$NEW3KTOT
# Load Benthic Data, for cells which are not developable for fish aqua set to NA
B.V <- rep(NA,times = length(F.V))
B.V[F.V > 0] <- sector_data.df$TOC.flux[Aqua.Full.Domain.Logical][F.V > 0]

# Run the eigenvector centrality diseaase model in MATLAB and then load the results.
# Write a .mat file with the filtered connectivity matrix
filename <- paste0(inpdatadir,"tmp.mat")
writeMat(filename,
eig = readMat(paste0(inpdatadir,
  'disease_connect_matrix.mat'))$disease.connect.matrix[F$Annuity > 0,F$Annuity > 0])
# Character vector to send to MATLAB from R. The function eigencentrality is derived from http://strategic.mit.edu/downloads.php?page=matlab_networks
code <- c("cd(strcat(pwd,\'/MSP_Model/Scripts/\'));",paste0('load \'',filename,'\';'),'d = abs(eigencentrality(eig));',
'save(\'tmp.mat\',\'d\')')
# Send arguments to matlab
run_matlab_code(code)
# Read the Mat file and remove the temporary one
D.V <- rep(NA,times = length(F.V))
D.V[F.V > 0] <- readMat(paste0(scrpdir,'tmp.mat'))$d
system2('rm',args = paste0(filename))
# Save all of the raw outputs of each sector model in a seperate file --> do later
print('Raw Impacts/Value.....')
Raw_Impacts <- data.frame(Mussel = M.V, Finfish = F.V, Kelp = K.V, Halibut = H.V,
  Viewshed_Mussel_Kelp = V_MK.V, Viewshed_Finfish = V_F.V, Benthic_Impacts = B.V,
  Disease_Risk = D.V) %>% glimpse()
# Make .mat files of the sector files
writeMat(paste0(inpdatadir,'Raw_Impacts.mat'),Raw_Impacts = Raw_Impacts)
writeMat(paste0(inpdatadir,'Raw_Impacts_FID.mat'),
      Raw_Impacts = select(mutate(Raw_Impacts,FID = fulldomain[Aqua.Full.Domain.Logical]),FID,Mussel,Finfish,Kelp,Halibut,Viewshed_Mussel_Kelp,Viewshed_Finfish,Benthic_Impacts,Disease_Risk))
writeMat(paste0(inpdatadir,'Raw_Impacts_FID.mat'),
      Raw_Impacts = select(mutate(Raw_Impacts,FID = fulldomain[Aqua.Full.Domain.Logical]),FID,Mussel,Finfish,Kelp,Halibut,Viewshed_Mussel_Kelp,Viewshed_Finfish,Benthic_Impacts,Disease_Risk))
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

X_n_i_p <- setNames(lapply(1:p, df = V_n_i_p, FUN = function(p,df){
      return(setNames(data.frame(t(do.call('rbind',lapply(1:length(true_sector_names), FUN = function(n){
        df[[p]][,n] / sum(apply(sapply(df,"[", ,n),MARGIN = 1, FUN = function(z){ifelse(!all(is.na(z)),max(z, na.rm = T),NA)}),na.rm = T)
      })))),true_sector_names))
  }),c('No_Development','Develop_Mussel','Develop_Finfish','Develop_Kelp'))
writeMat(file.path(paste0(outdatadir,'X_n_i_p.mat')),X_n_i_p = X_n_i_p)
# Create the sector weights (alpha's)
library(gtools)
epsilon <- .20 # Epsilon step size, default is 0.20
a_values <- seq(from = 0, to = 1, by = epsilon) # The unique values for each sector and site
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
  # # Save model results
  obj_i_mat <- as.matrix(obj_i)
  # Seperating the plans into more manageable sets
  obj_i_mat_split <- split(1:ncol(obj_i_mat), ceiling(seq_along(1:ncol(obj_i_mat))/93312))
  obj_i_mat.1 <- obj_i_mat[,obj_i_mat_split[[1]]]
  obj_i_mat.2 <- obj_i_mat[,obj_i_mat_split[[2]]]
  obj_i_mat.3 <- obj_i_mat[,obj_i_mat_split[[3]]]

  writeMat(file.path(paste0(outdatadir,'Policy_i_a_1.mat')), Policy_i_a_1 = as.matrix(obj_i_mat.1))
  writeMat(file.path(paste0(outdatadir,'Policy_i_a_2.mat')), Policy_i_a_2 = as.matrix(obj_i_mat.2))
  writeMat(file.path(paste0(outdatadir,'Policy_i_a_3.mat')), Policy_i_a_3 = as.matrix(obj_i_mat.3))

  write.table(as.matrix(obj_i_mat.1),file = file.path(paste0(outdatadir,'Policy_i_a_1.csv')), col.names = FALSE, row.names = FALSE)
  write.table(as.matrix(obj_i_mat.2),file = file.path(paste0(outdatadir,'Policy_i_a_2.csv')), col.names = FALSE, row.names = FALSE)
  write.table(as.matrix(obj_i_mat.3),file = file.path(paste0(outdatadir,'Policy_i_a_3.csv')), col.names = FALSE, row.names = FALSE)
  # system2('cp',args = c(inpdatadir,paste0(inpdatadir,'Lester_et_al_MSPsolutions_Evals_v3/')))
}else{
  print('loading planning results')
  # obj_i <- read.csv(file.path(paste0(outdatadir,'MSP_Planning_Results.csv')),header = FALSE)
  obj_i_load <- cbind(readMat(file.path(paste0(outdatadir,'Policy_i_a_1.mat')))[[1]],
                  readMat(file.path(paste0(outdatadir,'Policy_i_a_2.mat')))[[1]],
                    readMat(file.path(paste0(outdatadir,'Policy_i_a_3.mat')))[[1]])
  # readMat(file.path(paste0(outdatadir,'Policy_i_a_2.mat'))))
  # readMat(file.path(paste0(outdatadir,'Policy_i_a_3.mat')))
  # obj_i <- read.csv(file.path(paste0(outdatadir,'Policy_i_a.csv')),header = FALSE)
}
## Run the conventional model
# User defined function to scale each aqua sector by the impacted sector which is maximally impacted by the development of a specific aquaculture type
conventional.scaled <- function(z){
  tmp <- ifelse(!all(is.na(z)),max(z, na.rm = T),NA) # For a given site, find the maximumally impacted sector which is not NA
  return(ifelse(tmp != 0,tmp,1)) # Before returning, remove entries that are 0 in order to prevent from dividing by zero
}
# Apply the conventional.scaled function to each of the responses of the impacted sector for the three types of aquaculture
# and then divide by the individual aquaculture values
# C.df <- data.frame(M = Raw_Impacts$Mussel / apply(select(X_n_i_p[[1]],c(-1,-2,-3)) - select(X_n_i_p[[2]],c(-1,-2,-3)),MARGIN = 1, FUN = conventional.scaled),
#                 F = Raw_Impacts$Finfish / apply(select(X_n_i_p[[1]],c(-1,-2,-3)) - select(X_n_i_p[[3]],c(-1,-2,-3)), MARGIN = 1, FUN = conventional.scaled),
#                 K = Raw_Impacts$Kelp / apply(select(X_n_i_p[[1]],c(-1,-2,-3)) - select(X_n_i_p[[4]],c(-1,-2,-3)), MARGIN = 1, FUN = conventional.scaled)) %>% glimpse()
C.df <- data.frame(M = Raw_Impacts$Mussel / apply(select(X_n_i_p[[2]],c(-1,-2,-3)),MARGIN = 1, FUN = conventional.scaled),
                F = Raw_Impacts$Finfish / apply(select(X_n_i_p[[3]],c(-1,-2,-3)), MARGIN = 1, FUN = conventional.scaled),
                K = Raw_Impacts$Kelp / apply(select(X_n_i_p[[4]],c(-1,-2,-3)), MARGIN = 1, FUN = conventional.scaled)) %>% glimpse()
# Unconstrained conventional model
# Dataframe consisting of the FID, planning descision, and value of the aquaculture sector
U.C.Full.df.tmp <- cbind(1:i,fulldomain[Aqua.Full.Domain.Logical],t(apply(C.df, MARGIN = 1, FUN = function(x){
  x.ind <- which(x == max(x,na.rm = T))
  return(c(x.ind,x[x.ind]))
  })))
U.C.Full.df <- setNames(data.frame(U.C.Full.df.tmp[order(U.C.Full.df.tmp[,4],decreasing = T),]),
  c('i','FID','Policy','Value'))
# Using the U.C.Full.df dataframe creating a 1061 x 1062 dataframe of individual spatial plans
U.C.obj_i <- cbind(rep(0,i),sapply(1:i,FUN = function(x){
    tmp <- rep(0,i)
    tmp[U.C.Full.df$i[1:x]] <- U.C.Full.df$Policy[1:x]
    return(tmp)
  }))
write.csv(U.C.obj_i,file.path(paste0(inpdatadir,'U_C_obj_i.csv')),row.names = FALSE)
writeMat(file.path(paste0(inpdatadir,'U_C_obj_i.mat')),U_C_obj_i = U.C.obj_i)
# Constrained conventional model
if(readline(prompt = "Run Constrained Conventional Model? ") == 'Y'){
  C.df.tmp <- C.df %>%
              mutate(FID = fulldomain[Aqua.Full.Domain.Logical]) %>%
              select(FID, M, F, K) %>% glimpse()
  C.C.Full.df <- lapply(1:3,FUN = function(x){
    setNames(data.frame(cbind(1:i,fulldomain[Aqua.Full.Domain.Logical],C.df[,x])),c('i','FID','Value')) %>% mutate(Aqua_Type = c('Mussel','Finfish','Kelp')[x]) %>% arrange(desc(Value))
    # filter(setNames(data.frame(cbind(1:i,fulldomain[Aqua.Full.Domain.Logical],C.df[,x])),c('i','FID','Value')),Value > 0) %>% arrange(desc(Value))
    })
  C.C.obj_i <- cbind(rep(0,i),sapply(1:i,FUN = function(x){
    print(x)
    num_sites <- x # Number of sites to be developed in the spatial plan
    C.C.Full.df.tmp <- C.C.Full.df
    tmp <- rep(0,i) # empty spatial plan
    num_dev <- 0
    repeat{
      df <- filter(data.frame(rbindlist(lapply(C.C.Full.df.tmp,'[',1,))),Value > 0) # Make a data frame from the most suitable site of the three aquaculture forms, remove aqua sectors which are zero
      # print(nrow(df))
      current_choices <- c('Mussel','Finfish','Kelp')
      for(itor in 1:nrow(df)){
        df.dev <- df[which.max(df$Value),] # The first sector developed
        df.dev.FID <- df.dev$FID # The FID of the first aquaculture site to be developed in the triplet
        choice <- ifelse(df.dev$Aqua_Type != 'Mussel',ifelse(df.dev$Aqua_Type != 'Kelp',2,3),1)
        current_choices <- current_choices[which(current_choices != df.dev$Aqua_Type)]
        tmp[which(fulldomain[Aqua.Full.Domain.Logical] == df.dev.FID)] <- choice # Develop that site for the selected type of aquaculture
        C.C.Full.df.tmp <- lapply(C.C.Full.df.tmp,FUN = function(y){arrange(filter(y,FID != df.dev.FID),desc(Value))}) # Remove that FID from the ordered suitability for all three aqauculture types
        df <- filter(filter(data.frame(rbindlist(lapply(C.C.Full.df.tmp,'[',1,))),Value > 0),FID != df.dev.FID & Aqua_Type %in% current_choices) # Update the temporary dataframe for the current triplet
        if(length(which(tmp != 0)) == num_sites){break}
      }
      if(length(which(tmp != 0)) == num_sites){break}
    }
    spatial_plan <- tmp
    print(table(spatial_plan))
    return(spatial_plan)
    }))
  write.csv(C.C.obj_i,file.path(paste0(inpdatadir,'C_C_obj_i.csv')),row.names = FALSE)
  writeMat(file.path(paste0(inpdatadir,'C_C_obj_i.mat')),C_C_obj_i = C.C.obj_i)
}
# Run Halibut model for the MSP model and the two conventional model
# MSP
if(readline('Run Dynamic Halibut Model For MSP? Y/N ') == 'Y'){
  system2(matlab_root,
    args = c('-nodesktop','-noFigureWindows','-nosplash','-r',
    paste0("run\\(\\'",scrpdir,"Halibut/SCB_MSP_Dynamic.m\\'\\)")))
  system2(matlab_root,
    args = c('-nodesktop','-noFigureWindows','-nosplash','-r',
    paste0("run\\(\\'",scrpdir,"Halibut/MSP_SolutionPlots_v3_MSP.m\\'\\)")))
}
# Unconstrained
if(readline('Run Conventional Halibut For Conventional Models? ') == 'Y'){
  system2(matlab_root,
    args = c('-nodesktop','-noFigureWindows','-nosplash','-r',
    paste0("run\\(\\'",scrpdir,"Halibut/dynamic_eval_UC.m\\'\\)")))
  system2(matlab_root,
    args = c('-nodesktop','-noFigureWindows','-nosplash','-r',
    paste0("run\\(\\'",scrpdir,"Halibut/MSP_SolutionPlots_v3_UC.m\\'\\)")))
  # Constrained
  system2(matlab_root,
    args = c('-nodesktop','-noFigureWindows','-nosplash','-r',
    paste0("run\\(\\'",scrpdir,"Halibut/dynamic_eval_CC.m\\'\\)")))
  system2(matlab_root,
    args = c('-nodesktop','-noFigureWindows','-nosplash','-r',
    paste0("run\\(\\'",scrpdir,"Halibut/MSP_SolutionPlots_v3_CC.m\\'\\)")))
}
# Plot Data
U.C.obj_i <- read.csv(file.path(paste0(inpdatadir,'U_C_obj_i.csv')),stringsAsFactors = FALSE) #%>% select(-X)
C.C.obj_i <- read.csv(file.path(paste0(inpdatadir,'C_C_obj_i.csv')),stringsAsFactors = FALSE) #%>% select(-X)
true_sector_names <- c('Mussel','Finfish','Kelp','Halibut','Viewshed','Benthic','Disease')
Static.values.data <- setNames(data.frame(sapply(readMat('~/MSP_Model/Input/Data/EFPayoff_a_X_wrt_DM.mat'),FUN = t)),true_sector_names)
Unconstrained_data <- setNames(data.frame(sapply(readMat('~/MSP_Model/Output/Data/EFPayoff_a_X_wrt_DM_UC.mat'),FUN = function(x){t(ifelse(is.na(x),0,x))})) %>%
  select(EFPayoff.a.M.wrt.DM,EFPayoff.a.F.wrt.DM,EFPayoff.a.K.wrt.DM,EFPayoff.a.H.wrt.DM,EFPayoff.a.V.wrt.DM,EFPayoff.a.B.wrt.DM,EFPayoff.a.D.wrt.DM),true_sector_names)
Constrained_data <- setNames(data.frame(sapply(readMat('~/MSP_Model/Output/Data/EFPayoff_a_X_wrt_DM_CC.mat'),FUN = function(x){t(ifelse(is.na(x),0,x))})) %>%
  select(EFPayoff.a.M.wrt.DM,EFPayoff.a.F.wrt.DM,EFPayoff.a.K.wrt.DM,EFPayoff.a.H.wrt.DM,EFPayoff.a.V.wrt.DM,EFPayoff.a.B.wrt.DM,EFPayoff.a.D.wrt.DM),true_sector_names)
Master.matrix.max <- rbind(Static.values.data,Unconstrained_data,Constrained_data)
color.vector.max=c(rep('coral',length.out=nrow(Static.values.data)),rep('purple',length.out=nrow(Unconstrained_data)),rep('green',length.out=nrow(Constrained_data)))


# Static.values.data %>% select(Mussel,Finfish,Kelp,Viewshed) %>%
#   filter(Viewshed > .95 & Mussel > 0 & Finfish > 0 & Kelp > 0) %>%
#   distinct() %>% mutate(Cumulative = apply(Static.values.data %>%
#   select(Mussel,Finfish,Kelp,Viewshed) %>% filter(Viewshed > .95 & Mussel > 0 & Finfish > 0 & Kelp > 0) %>%
#   distinct() %>% select(Mussel, Finfish, Kelp),MARGIN = 1,FUN = mean)) %>% arrange(desc(Cumulative))

# Static.values.data %>% select(Finfish,Disease) %>% filter(Disease > .99) %>% arrange(desc(Finfish)) %>% distinct()

seeds <- setNames(data.frame(t(readMat('~/MSP_Model/Input/Data/EFPayoff_a_X_wrt_DM_filter.mat')[[1]])),true_sector_names) %>%
  select(Mussel, Finfish, Kelp) %>%
  mutate(Mussel = Mussel * sum(Raw_Impacts$Mussel)) %>%
  mutate(Finfish = Finfish * sum(Raw_Impacts$Finfish)) %>%
  mutate(Kelp = Kelp * sum(Raw_Impacts$Kelp)) %>% distinct() %>% summarise_all(funs(min,max)) %>% select(Mussel_min, Mussel_max, Finfish_min, Finfish_max, Kelp_min, Kelp_max)

Static.values.data %>% select(Mussel, Finfish, Kelp) %>% mutate(Cumliative = Mussel + Finfish + Kelp) %>% distinct() %>% filter(Cumliative == max(Cumliative))
width=7
height=5.5
res=1000
units='in'
text.size<-12
patch.size=1.5
cols <- c("Mussel"="Blue","Finfish"="Salmon","Kelp"="Green","Halibut"="Burlywood","Viewshed"='Cyan',"Benthic"='Orange',"Disease"='Black')
formatList <- c('png','eps','pdf') # Specify formats to plot

# Custom plotting functions (used throughout figure code)
 png_load <- function(imageList, theme, labs, subtitleList = NULL, barsObj = NULL){
    images <- list()
    if(is.null(subtitleList)){
        subtitleList <- c('A','B','C','D')
    }
    for(image_itor in 1:length(imageList)){
        img <- readPNG(paste0(inpfigdir,imageList[image_itor]),native=T,info=T)
        g <- rasterGrob(img, interpolate=TRUE)
        ggObj<-qplot(1:10, 1:10, geom="blank") +
            annotation_custom(g, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +
            ggtitle(subtitleList[image_itor]) +
            theme + labs
        images[[image_itor]] <- ggObj
    }
    return(images)
}

figure_output <- function(formatList, figure_object, outfigdir, file_name, width, height, units){
    for(itor in 1:length(formatList)){
        dir.create(paste0(outfigdir, 'Main_MS/'), showWarnings = FALSE)
        folder <- paste0(outfigdir, 'Main_MS/', formatList[itor])
        dir.create(folder, showWarnings = FALSE)
        file_name_full <- file.path(folder, paste0(file_name, '.', formatList[itor]))
        print(file_name_full)
        ggsave(filename = file_name_full,
            plot = figure_object,
            width = width,
            height = height,
            units = units,
            dpi = 1000)
    }
    graphics.off()
}
# Theme used throughout primary text unless otherwise specified.
theme = theme(plot.margin = unit(c(.2,0,.2,1), units = "lines"),
                axis.text = element_blank(),
                axis.title = element_blank(),
                axis.ticks = element_blank(),
                axis.ticks.length = unit(0, "lines"),
                axis.ticks.margin = unit(0, "lines"),
                panel.background=element_rect(fill="white"),
                panel.grid=element_blank(),
                plot.title=element_text(hjust=0))
labs = labs(x = NULL, y = NULL)

# Figure 1
imageList=c("Fig1A Capture.png",
    "MusselValueApril.png",
    "FishValueApril.png",
    "KelpValueApril.png")
fig1List = png_load(imageList,theme,labs)
fig1 <- arrangeGrob(grobs = fig1List,
    ncol=2,
    nrow=2,
    padding = unit(-.5,'lines'))
figure_output(formatList, fig1, outfigdir, file_name='Fig 1', width, height, units)
# Figure 2
source('~/MSP_Model/Scripts/Fig2.r')
figure2(formatList)
# dev.off()
# Figure 3
imageList=c("Figure_3_ALL.png",
    "Figure_3_mussel.png",
    "Figure_3_finfish.png",
    "Figure_3_kelp.png")
fig3List = png_load(imageList,theme,labs)
fig3 <- arrangeGrob(grobs = fig3List,
    ncol=2,
    nrow=2,
    padding = unit(-.5,'lines'))
figure_output(formatList, fig3, outfigdir, file_name='Fig 3', width, height, units)

# Figure 4
EFPayoff_filter <- t(readMat('~/MSP_Model/Input/Data/EFPayoff_a_X_wrt_DM_filter.mat')[[1]])
Low.impact.solutions=setNames(data.frame(EFPayoff_filter),c('Mussel','Finfish','Kelp','Halibut','Viewshed','Benthic','Disease'))
LI.names=names(Low.impact.solutions)
ID=seq(to=nrow(Low.impact.solutions),from=1,by=1)
size.tmp=dim(Low.impact.solutions)
Mussel.LI=Low.impact.solutions[,1]
Finfish.LI=Low.impact.solutions[,2]
Kelp.LI=Low.impact.solutions[,3]
Halibut.LI=Low.impact.solutions[,4]
View.LI=Low.impact.solutions[,5]
Benthic.LI=Low.impact.solutions[,6]
Disease.LI=Low.impact.solutions[,7]
Sector.LI=c(rep(LI.names[1],times=size.tmp[1]),
            rep(LI.names[2],times=size.tmp[1]),
            rep(LI.names[3],times=size.tmp[1]),
            rep(LI.names[4],times=size.tmp[1]),
            rep(LI.names[5],times=size.tmp[1]),
            rep(LI.names[6],times=size.tmp[1]),
            rep(LI.names[7],times=size.tmp[1]))
for(itor in 1:ncol(Low.impact.solutions)){
  ID.LI.tmp=1:size.tmp[1]
  if(itor>1){
    ID.LI=c(ID.LI,ID.LI.tmp)
  }else{
    ID.LI=ID.LI.tmp
  }
}
value.LI.tmp=c(Mussel.LI,Finfish.LI,Kelp.LI,Halibut.LI,View.LI,Benthic.LI,Disease.LI)
Low.impact.solutions=data.frame(ID.LI,value.LI.tmp,Sector.LI)
names(Low.impact.solutions)=c('ID','Value','Sector')
Low.impact.solutions$Sector=factor(Low.impact.solutions$Sector, levels=c('Mussel','Finfish','Kelp','Halibut','Viewshed','Benthic','Disease'))
p.bar<-ggplot(data = Low.impact.solutions,aes(x=ID,y=Value,fill=Sector,color=Sector))+
  geom_bar(stat="identity")+facet_grid(.~Sector)+ggtitle('A')+
  ylab('Value [% of maximum]') +
  scale_y_continuous(labels = percent,limits=c(0,1))+
  scale_fill_manual(values=cols)+
  scale_color_manual(values=cols)+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        panel.background=element_rect(color='white',fill='white'),panel.spacing = unit(0, "lines"),
        strip.background=element_rect(fill='white'),strip.text=element_text(size=text.size*.75),
        panel.grid=element_blank(),axis.title.y=element_text(size=text.size,color="black"),
        axis.text.y=element_text(size=text.size,color="black"),legend.title=element_text(size=text.size*.75),
        legend.text=element_text(size=.75*text.size),plot.title=element_text(hjust=0),plot.margin = unit(c(.2,.2,.2,.2), units = "lines"))

img_case_study_percent <- readPNG(paste0(inpfigdir,'Figure4B.png'),native=T,info=T)
g_case_study_percent <- rasterGrob(img_case_study_percent, interpolate=TRUE,just='center')

p.map<-qplot(1:10, 1:10, geom="blank") + ggtitle('B') +
  annotation_custom(g_case_study_percent, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)+
  theme(plot.margin = unit(c(.2,.2,.2,.2), units = "lines"),axis.text.x = element_blank(),axis.text.y = element_blank(),axis.ticks = element_blank(),
        axis.title=element_blank(),panel.background=element_rect(fill="white"),panel.grid=element_blank(),
        plot.title=element_text(hjust=0)) + theme

fig4 <- arrangeGrob(p.bar,p.map,ncol=2,nrow=1,padding=unit(-.5,'line'))
figure_output(formatList, fig4, outfigdir, file_name='Fig 4', width=width * 2, height=height,units = units)

# ggsave(paste0(outfigdir,'Fig 4.png'),fig4,width=width * 2, height=height,units = units)

# Figure 5
print('foo')
EFPayoff_a_X_wrt_DM_SeedS1S2S3 <- setNames(data.frame(t(readMat(paste0(inpdatadir,'EFPayoff_a_X_wrt_DM_Seed_bc_set.mat'))[[1]])),
        c('Mussel','Finfish','Kelp','Halibut','Viewshed','Benthic','Disease')) %>% mutate(Seed = c('S1','S2','S3','S4','S5')) %>% gather(Sector,Value,-Seed)
p.bar <- ggplot(EFPayoff_a_X_wrt_DM_SeedS1S2S3,aes(x=factor(Sector, levels = c('Mussel','Finfish','Kelp','Halibut','Viewshed','Benthic','Disease')),y=Value,fill=factor(Seed)))+
  geom_bar(stat="identity",position="dodge" )+
  ylab('Value [% of maximum]') +
  scale_fill_discrete(name="Seed",labels=c('Prioritize Existing \nSectors','Prioritize Aquaculture \nSectors','Balance Existing \nSectors and Aquaculture')) +
  ggtitle('A') + scale_y_continuous(labels = scales::percent) +
  theme(axis.title.x=element_blank(),
        axis.ticks.x=element_blank(),
        panel.background=element_rect(color='white',fill='white'),
        strip.background=element_rect(fill='white'),strip.text=element_text(size=text.size),
        panel.grid=element_blank(),axis.title.y=element_text(size=text.size,color="black"),
        axis.text.x=element_text(size=text.size*.65,color="black"),axis.text.y=element_text(size=text.size*.65,color="black"),
        legend.title=element_blank(),
        legend.text=element_text(size=text.size * .45),
        plot.title=element_text(hjust=0),legend.position = 'bottom',
        legend.key.size = unit(1,"lines"),
        legend.key = element_rect(size = 3))
imageList <- c("Seed_Plan1.png","Seed_Plan3.png","Seed_Plan5.png")
fig5List <- png_load(imageList,theme,labs,subtitleList=c('B','C','D'))
# fig5list <- c(p.bar, fig5List)
# print(fig5List)
# fig5List[1] <- p.bar
# fig5List <- c(p.bar, fig5List[[1]], fig5List[2], fig5List[3])
fig5 <- arrangeGrob(p.bar, fig5List[[1]], fig5List[[2]], fig5List[[3]],ncol=2,nrow=2,layout_matrix = rbind(c(1,2),c(3,4)),padding=unit(-1,'line'))
figure_output(formatList, fig5, outfigdir, file_name='Fig 5', width=width*1.15, height=height*1.15, units)
# ggsave(paste0(outfigdir,'Fig 5.png'),fig5,width=width*1.15, height=height*1.15,units = units)
# Figure 6
    source(paste0(scrpdir,'value.of.MSP.loop_interpol.R'))
    source(paste0(scrpdir,'value.of.MSP.fx_2_interpol.R'))
    source(paste0(scrpdir,'figure_5_code.R'))
    MSP.value.data=rbind(Mussel.value.tmp,Finfish.value.tmp,Kelp.value.tmp)
    MSP.value.data.points=rbind(Mussel.value.tmp.points,Finfish.value.tmp.points,Kelp.value.tmp.points)
    names(MSP.value.data)=c('Aquaculture.Value','Value.of.MSP','Group','Sector.Name','Sector.Type','Type.of.Conventional')
    names(MSP.value.data.points)=c('Aquaculture.Value','Value.of.MSP','Group','Sector.Name','Sector.Type','Type.of.Conventional')
    MSP.value.data$Sector.Name=factor(MSP.value.data$Sector.Name,levels = c('Mussel','Finfish','Kelp','Halibut','Viewshed','Benthic','Disease'))
    Value.of.MSP.grid.plot<-ggplot(data=MSP.value.data)+
      geom_point(data=subset(MSP.value.data[as.integer(MSP.value.data$Group)!=as.integer(MSP.value.data$Sector.Name),]),
                 aes(x=Aquaculture.Value,y=Value.of.MSP,shape=Type.of.Conventional,color=Type.of.Conventional),size=1.25)+
      facet_grid(Sector.Name~Group)+scale_y_continuous(labels = percent,limits=c(0,1),breaks=c(0,.5,1))+
      scale_x_continuous(labels = percent,breaks=c(0,.5,1))+
      geom_text(data=subset(MSP.value.data[as.integer(MSP.value.data$Group)==as.integer(MSP.value.data$Sector.Name),]),x=.5,y=.5,size=15,label='NA')+
      geom_point(data=MSP.value.data.points,aes(x=Aquaculture.Value,y=Value.of.MSP,shape=Type.of.Conventional,color=Type.of.Conventional),size=1.2)+
      xlab("Aquaculture Value")+ylab("Value of Marine Spatial Planning")+
      scale_colour_manual(name = "Conventional Planning :",labels=c("Constrained","Unconstrained"),values=c("coral1","mediumorchid1"))+
      scale_shape_manual(name = "Conventional Planning :",labels=c("Constrained","Unconstrained"),values=c(16,24))+
      scale_fill_manual(name = "Conventional Planning :",labels=c("Constrained","Unconstrained"),values=c("coral1","mediumorchid1"))+
      theme_bw(base_size = 15)+theme(panel.grid=element_blank(),legend.position="bottom")

fig6 <- Value.of.MSP.grid.plot+theme(axis.title=element_text(size=12),panel.spacing = unit(1, "lines"),
                             strip.background=element_rect(fill='white'),strip.text=element_text(size=10),axis.ticks.margin=unit(1,'lines'))
figure_output(formatList, fig6, outfigdir, file_name='Fig 6', width=8, height=8, units)
# ggsave(paste0(outfigdir,'Fig 6.png'),fig6,width=8, height=8,units=units)
## Figure 1
## SI Figures
# current.directory.code <- '~/Desktop/Code/SI Figures Code/'
# current.directory.figures <- '~/Desktop/Code/SI Figures/'
# Coastline.data<-read.csv(file='Coastline.csv',header=F)
# plot.margin = unit(c(.50,0,-.15,1), units = "lines")
theme_3 = theme(plot.margin = unit(c(.1,.1,-.25,.1), units = "lines"),
              axis.text = element_blank(),
              axis.title = element_blank(),
              axis.ticks = element_blank(),
              # axis.ticks.length = unit(0, "lines"),
              # axis.ticks.margin = unit(0, "lines"),
              panel.background=element_rect(fill="white"),
              panel.grid=element_blank(),
              plot.title=element_text(hjust =.37,vjust=0,margin=margin(0,0,0,0)))

theme_3A = theme(plot.margin = unit(c(.1,.1,-.25,.1), units = "lines"),
              axis.text = element_blank(),
              axis.title = element_blank(),
              axis.ticks = element_blank(),
              # axis.ticks.length = unit(0, "lines"),
              # axis.ticks.margin = unit(0, "lines"),
              panel.background=element_rect(fill="white"),
              panel.grid=element_blank(),
              plot.title=element_text(hjust =.37,vjust=0,margin=margin(0,0,0,0)))

theme_3A = theme(plot.margin = unit(c(.1,.1,-.25,.1), units = "lines"),
              axis.text = element_blank(),
              axis.title = element_blank(),
              axis.ticks = element_blank(),
              # axis.ticks.length = unit(0, "lines"),
              # axis.ticks.margin = unit(0, "lines"),
              panel.background=element_rect(fill="white"),
              panel.grid=element_blank(),
              plot.title=element_text(hjust =.37,vjust=0,margin=margin(0,0,0,0)))

theme_3B = theme(plot.margin = unit(c(-.25,.1,-.25,.1), units = "lines"),
              axis.text = element_blank(),
              axis.title = element_blank(),
              axis.ticks = element_blank(),
              # axis.ticks.length = unit(0, "lines"),
              # axis.ticks.margin = unit(0, "lines"),
              panel.background=element_rect(fill="white"),
              panel.grid=element_blank(),
              plot.title=element_text(hjust =.37,vjust=0,margin=margin(0,0,0,0)))

theme_3C = theme(plot.margin = unit(c(-.25,.1,.1,.1), units = "lines"),
              axis.text = element_blank(),
              axis.title = element_blank(),
              axis.ticks = element_blank(),
              # axis.ticks.length = unit(0, "lines"),
              # axis.ticks.margin = unit(0, "lines"),
              panel.background=element_rect(fill="white"),
              panel.grid=element_blank(),
              plot.title=element_text(hjust =.37,vjust=0,margin=margin(0,0,0,0)))

theme_2A = theme(plot.margin = unit(c(.1,.1,-.25,.1), units = "lines"),
              axis.text = element_blank(),
              axis.title = element_blank(),
              axis.ticks = element_blank(),
              # axis.ticks.length = unit(0, "lines"),
              # axis.ticks.margin = unit(0, "lines"),
              panel.background=element_rect(fill="white"),
              panel.grid=element_blank(),
              plot.title=element_text(hjust =.27,vjust=0,margin=margin(0,0,0,0)))
theme_2B = theme(plot.margin = unit(c(-.25,.1,.1,.1), units = "lines"),
              axis.text = element_blank(),
              axis.title = element_blank(),
              axis.ticks = element_blank(),
              # axis.ticks.length = unit(0, "lines"),
              # axis.ticks.margin = unit(0, "lines"),
              panel.background=element_rect(fill="white"),
              panel.grid=element_blank(),
              plot.title=element_text(hjust =.27,vjust=0,margin=margin(0,0,0,0)))

labs = labs(x = NULL, y = NULL)
si_figure_output <- function(formatList, figure_object, outfigdir, file_name, width, height, units){
    for(itor in 1:length(formatList)){
        dir.create(paste0(outfigdir, 'SI/'), showWarnings = FALSE)
        folder <- paste0(outfigdir, 'SI/', formatList[itor])
        dir.create(folder, showWarnings = FALSE)
        file_name_full <- file.path(folder, paste0(file_name, '.', formatList[itor]))
        print(file_name_full)
        ggsave(filename = file_name_full,
            plot = figure_object,
            width = width,
            height = height,
            units = units,
            dpi = 1000)
    }
    graphics.off()
}
# Figure S1
img_s1 <- readPNG(paste0(inpfigdir,'Fig S1.png'),native=T,info=T)
g_s1 <- rasterGrob(img_s1, interpolate=TRUE)

S1 <- qplot(1:10, 1:10, geom="blank") +
      annotation_custom(g_s1, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)+
      theme_3 + labs
si_figure_output(formatList, S1, outfigdir, file_name='Fig S1', width=width, height=height,units = units)
# ggsave(paste0(outfigdir,'Fig S1.png'),S1,width=width, height=height,units = units)
# Figure S2
img_S2A <- readPNG(paste0(inpfigdir,'S2A.png'),native=T,info=T)
g_S2A <- rasterGrob(img_S2A, interpolate=TRUE)

S2A<-qplot(1:10, 1:10, geom="blank") + ggtitle('A') +
      annotation_custom(g_S2A, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)+
      theme_3A + labs #gin = unit(c(0,0,0,0),'lines'))

img_S2B <- readPNG(paste0(inpfigdir,'S2B.png'),native=T,info=T)
g_S2B <- rasterGrob(img_S2B, interpolate=TRUE)

S2B<-qplot(1:10, 1:10, geom="blank") + ggtitle('B') +
      annotation_custom(g_S2B, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)+
      theme_3B + labs #+ theme(plot.margin = unit(c(0,0,0,0),'lines'))

img_S2C <- readPNG(paste0(inpfigdir,'S2C.png'),native=T,info=T)
g_S2C <- rasterGrob(img_S2C, interpolate=TRUE)

S2C<-qplot(1:10, 1:10, geom="blank") + ggtitle('C') +
  annotation_custom(g_S2C, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)+
  theme_3C + labs #+ theme(plot.margin = unit(c(0,0,0,0),'lines'))

S2 <- arrangeGrob(S2A,S2B,S2C,ncol=1,nrow=3,padding=unit(0,'line'))
si_figure_output(formatList, S2, outfigdir, file_name='Fig S2', width=width, height=height,units = units)

# ggsave(paste0(outfigdir,'Fig S2.png'),S2, width = width, height=height,units=units)
# S3
img_S3A <- readPNG(paste0(inpfigdir,'S3A.png'),native=T,info=T)
g_S3A <- rasterGrob(img_S3A, interpolate=TRUE)

S3A<-qplot(1:10, 1:10, geom="blank") + ggtitle('A') +
  annotation_custom(g_S3A, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)+
  theme_3A + labs #+ #theme(plot.margin = unit(c(1,0,-.15,1),'lines'))

img_S3B <- readPNG(paste0(inpfigdir,'S3B.png'),native=T,info=T)
g_S3B <- rasterGrob(img_S3B, interpolate=TRUE)

S3B<-qplot(1:10, 1:10, geom="blank") + ggtitle('B') +
  annotation_custom(g_S3B, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)+
  theme_3B + labs

img_S3C <- readPNG(paste0(inpfigdir,'S3C.png'),native=T,info=T)
g_S3C <- rasterGrob(img_S3C, interpolate=TRUE)

S3C<-qplot(1:10, 1:10, geom="blank") + ggtitle('C') +
  annotation_custom(g_S3C, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)+
  theme_3C + labs
S3 <- arrangeGrob(S3A,S3B,S3C,ncol=1,nrow=3,padding=unit(0,'line'))
si_figure_output(formatList, S3, outfigdir, file_name='Fig S3', width=width, height=height,units = units)
# ggsave(paste0(outfigdir,'Fig S3.png'),S3,width=width, height=height,units=units)
# S4
img_S4A <- readPNG(paste0(inpfigdir,'S5A.png'),native=T,info=T)
g_S4A <- rasterGrob(img_S4A, interpolate=TRUE)

S4A<-qplot(1:10, 1:10, geom="blank") + ggtitle('A') +
  annotation_custom(g_S4A, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)+
  theme_2A + labs #+ theme(plot.margin = unit(c(1,0,-.15,1),'lines'))

img_S4B <- readPNG(paste0(inpfigdir,'S5B.png'),native=T,info=T)
g_S4B <- rasterGrob(img_S4B, interpolate=TRUE)

S4B<-qplot(1:10, 1:10, geom="blank") + ggtitle('B') +
  annotation_custom(g_S4B, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)+
  theme_2B + labs

S4 <- arrangeGrob(S4A,S4B,ncol=1,nrow=2,padding=unit(-.5,'line'))
si_figure_output(formatList, S4, outfigdir, file_name='Fig S4', width=width, height=height,units = units)

ggsave(paste0(outfigdir,'Fig S4.png'),S4,width=width, height=height,units=units)
# S5
img_S5A <- readPNG(paste0(inpfigdir,'S6A.png'),native=T,info=T)
g_S5A <- rasterGrob(img_S5A, interpolate=TRUE)

S5A<-qplot(1:10, 1:10, geom="blank") + ggtitle('A') +
  annotation_custom(g_S5A, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)+
  theme_2A + labs #+ theme(plot.margin = unit(c(1,0,-.15,1),'lines'))

img_S5B <- readPNG(paste0(inpfigdir,'S6B.png'),native=T,info=T)
g_S5B <- rasterGrob(img_S5B, interpolate=TRUE)

S5B<-qplot(1:10, 1:10, geom="blank") + ggtitle('B') +
  annotation_custom(g_S5B, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)+
  theme_2B + labs

S5 <- arrangeGrob(S5A,S5B,ncol=1,nrow=2,padding=unit(-.5,'line'))
ggsave(paste0(outfigdir,'Fig S5.png'),S5,width=width, height=height,units=units)
# S6
img_S6A <- readPNG(paste0(inpfigdir,'S7A.png'),native=T,info=T)
g_S6A <- rasterGrob(img_S6A, interpolate=TRUE)

S6A<-qplot(1:10, 1:10, geom="blank") + ggtitle('A') +
  annotation_custom(g_S6A, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)+
  theme(axis.text = element_blank(),
                axis.title = element_blank(),
                axis.ticks = element_blank(),
                # axis.ticks.length = unit(0, "lines"),
                # axis.ticks.margin = unit(0, "lines"),
                panel.background=element_rect(fill="white"),
                panel.grid=element_blank(),
                plot.title=element_text(hjust=.1,vjust=0,margin=margin(0,0,0,0)),
                plot.margin = unit(c(.1,.1,-.25,-.25), units = "lines")) + labs

img_S6B <- readPNG(paste0(inpfigdir,'S7B.png'),native=T,info=T)
g_S6B <- rasterGrob(img_S6B, interpolate=TRUE)

S6B<-qplot(1:10, 1:10, geom="blank") + ggtitle('B') +
  annotation_custom(g_S6B, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)+
  theme(axis.text = element_blank(),
                axis.title = element_blank(),
                axis.ticks = element_blank(),
                # axis.ticks.length = unit(0, "lines"),
                # axis.ticks.margin = unit(0, "lines"),
                panel.background=element_rect(fill="white"),
                panel.grid=element_blank(),
                plot.title=element_text(hjust=.1,vjust=0,margin=margin(0,0,0,0)),
                plot.margin = unit(c(.1,-.25,.1,-.25), units = "lines")) + labs
# plot.margin = unit(c(0,0,0,0), units = "lines"),
img_S6C <- readPNG(paste0(inpfigdir,'S7C.png'),native=T,info=T)
g_S6C <- rasterGrob(img_S6C, interpolate=TRUE)

S6C<-qplot(1:10, 1:10, geom="blank") + ggtitle('C') +
  annotation_custom(g_S6C, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)+
  theme(axis.text = element_blank(),
                axis.title = element_blank(),
                axis.ticks = element_blank(),
                # axis.ticks.length = unit(0, "lines"),
                # axis.ticks.margin = unit(0, "lines"),
                panel.background=element_rect(fill="white"),
                panel.grid=element_blank(),
                plot.title=element_text(hjust=.1,vjust=0,margin=margin(0,0,0,0)),
                plot.margin = unit(c(-.25,.1,.1,-.25), units = "lines")) + labs
# plot.margin = unit(c(0,0,0,0), units = "lines"),
img_S6D <- readPNG(paste0(inpfigdir,'S7D.png'),native=T,info=T)
g_S6D <- rasterGrob(img_S6D, interpolate=TRUE)

S6D<-qplot(1:10, 1:10, geom="blank") + ggtitle('D') +
  annotation_custom(g_S6D, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)+
  theme(axis.text = element_blank(),
                axis.title = element_blank(),
                axis.ticks = element_blank(),
                # axis.ticks.length = unit(0, "lines"),
                # axis.ticks.margin = unit(0, "lines"),
                panel.background=element_rect(fill="white"),
                panel.grid=element_blank(),
                plot.title=element_text(hjust=.1,vjust=0,margin=margin(0,0,0,0)),
                plot.margin = unit(c(-.25,-.25,.1,.1), units = "lines")) + labs
# plot.margin = unit(c(0,0,0,0), units = "lines"),
S6 <- arrangeGrob(S6A,S6B,S6C,S6D,ncol=2,nrow=2,padding=unit(-.5,'line'))
ggsave(paste0(outfigdir,'Fig S6.png'),S6,width=width, height=height,units=units)

# S7
U.C.Summary <- rbindlist(lapply(1:ncol(U.C.obj_i),FUN = function(x){
                  data.frame(Mussel = length(which(U.C.obj_i[,x] == 1)),Finfish = length(which(U.C.obj_i[,x] == 2)),Kelp = length(which(U.C.obj_i[,x] == 3)))
                })) %>% mutate(Itor = 1:ncol(U.C.obj_i)) %>% gather(Sector, Number, -Itor) %>% glimpse()
C.C.Summary <- rbindlist(lapply(1:ncol(C.C.obj_i),FUN = function(x){
                  data.frame(Mussel = length(which(C.C.obj_i[,x] == 1)),Finfish = length(which(C.C.obj_i[,x] == 2)),Kelp = length(which(C.C.obj_i[,x] == 3)))
                })) %>% mutate(Itor = 1:ncol(U.C.obj_i)) %>% gather(Sector, Number, -Itor) %>% glimpse()
S7A <- ggplot() + geom_line(data = U.C.Summary, aes(x = Itor, y = Number, group = factor(Sector), color = factor(Sector))) +
        scale_fill_discrete(name="Seed",labels=c('Mussel','Finfish','Kelp')) +
        ggtitle('A') +
        xlab('Total Number of Sites Developed') +
        ylab('Number of Sites Developed') +
        theme(panel.background=element_rect(color='white',fill='white'),panel.spacing = unit(.75, "lines"),
              strip.background=element_rect(fill='white'),strip.text=element_text(size=text.size*.75),
              panel.grid=element_blank(),axis.title.y=element_text(size=text.size,color="black"),
              axis.text.y=element_text(size=text.size,color="black"),legend.title=element_blank(),
              axis.title.x=element_text(size=text.size,color="black"),
              axis.text.x=element_text(size=text.size,color="black"),
              legend.text=element_text(size=text.size),plot.title=element_text(hjust=0),
              legend.position = c(.17, .80),
              panel.border = element_rect(colour = "black", fill=NA, size=1),legend.key = element_rect(color='white',fill='white'))
S7B <- ggplot() + geom_line(data = C.C.Summary, aes(x = Itor, y = Number, group = factor(Sector), color = factor(Sector))) +
        scale_fill_discrete(name="Seed",labels=c('Mussel','Finfish','Kelp')) +
        ggtitle('B') +
        xlab('Total Number of Sites Developed') +
        ylab('Number of Sites Developed') +
        theme(panel.background=element_rect(color='white',fill='white'),panel.spacing = unit(.75, "lines"),
              strip.background=element_rect(fill='white'),strip.text=element_text(size=text.size*.75),
              panel.grid=element_blank(),axis.title.y=element_text(size=text.size,color="black"),
              axis.text.y=element_text(size=text.size,color="black"),legend.title=element_blank(),
              axis.title.x=element_text(size=text.size,color="black"),
              axis.text.x=element_text(size=text.size,color="black"),
              legend.text=element_text(size=text.size),plot.title=element_text(hjust=0),
              legend.position = c(.17, .80),
              panel.border = element_rect(colour = "black", fill=NA, size=1),legend.key = element_rect(color='white',fill='white'))
S7 <- arrangeGrob(S7A, S7B, ncol = 1, nrow = 2)
ggsave(paste0(outfigdir,'Fig S7.png'),S7,width=7, height=7,units=units)
