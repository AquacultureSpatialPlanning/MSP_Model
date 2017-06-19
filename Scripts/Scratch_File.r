# setwd('MSP_Model/')
# source('R_Libraries.r')
# choice <- ifelse(readline(prompt = 'Install? y/n') == 'Y',TRUE,FALSE)
# print(choice)
# R_Libraries(choice) # After the first initial run this can be set to F
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
# Remove unprofitable sites and generate seperate vectors for
# Those sites in which ventures will be profitable for at least one type
# of aquaculture (var Aqua.Full.Domain),
Aqua.Full.Domain <- data.frame(M$Annuity,F$Annuity,K$Annuity)
Aqua.Full.Domain.Logical <- apply(1 * (Aqua.Full.Domain > 0),FUN=sum,MARGIN=1) > 0
# Figure 3
Fig3 <- setNames(data.frame(do.call('cbind',readMat('~/MSP_Model/Input/Data/Lester_et_al_MSPsolutions_Evals_v3/Fig3_data.mat'))),c('Aquaculture','Finfish','Kelp','Mussel')) %>%
        mutate(FID = fulldomain[Aqua.Full.Domain.Logical]) %>%
        select(FID, Aquaculture, Mussel, Finfish, Kelp) %>%
        write.csv('~/MSP_Model/Figure_3_MS.csv',row.names = FALSE)
# Figure 4
Fig4 <- setNames(data.frame(FID = fulldomain[Aqua.Full.Domain.Logical], Freq = readMat('~/MSP_Model/Input/Data/Lester_et_al_MSPsolutions_Evals_v3/Fig4_data.mat')[[1]]),c('FID','Freq')) %>%
        write.csv('~/MSP_Model/Figure_4_MS.csv',row.names = FALSE)
# Figure 5
Fig5 <- setNames(data.frame(do.call('cbind',readMat('~/MSP_Model/Input/Data/Lester_et_al_MSPsolutions_Evals_v3/Policy_i_a_SX.mat')) - 1),c('S1','S2','S3')) %>%
        mutate(FID = fulldomain[Aqua.Full.Domain.Logical]) %>%
        select(FID, S1, S2, S3) %>%
        write.csv('~/MSP_Model/Figure_5_MS.csv')
