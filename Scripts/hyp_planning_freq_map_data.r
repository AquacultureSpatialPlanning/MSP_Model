current.directory.data='~/Desktop/Code/MSP Planning Results April 2016'
# current.directory.mac='/Users/joelstevens/Desktop/December 13th work/Code/MSP Planning Results April 2016'
wkdir <- getwd()
# install.packages(c('TeachingDemos','maps','mapdata','maptools','scales','ggmap','ggplot2','grid','GGally','gridExtra','jpeg','R.matlab','png','shape','DDHFm','tiff'))
source('~/Useful_R_packages.r')
Useful_R_packages(F,F)
require(colorout)
library(TeachingDemos)
library(maps)
library(mapdata)
library(maptools)
library(scales)
library(ggmap)
library(ggplot2)
library(grid)
library(GGally)
library(gridExtra)
library(jpeg)
library(R.matlab)
library(png)
library(shape)
library(DDHFm)
library(tiff)
library(dplyr)

values <- read.csv(file='~/Desktop/Code/MSP Planning Results April 2016/Dynamic_Values_Export.csv', header = FALSE)
colnames(values)=c('Mussel','Finfish','Kelp','Halibut','Viewshed','Benthic','Disease')
I <- which(values$Mussel>=.05&values$Finfish>=.05&values$Kelp>=.05&
                           values$Halibut>=.95&values$Viewshed>=.95&values$Benthic>=.95&
                           values$Disease>=.95)
# Save a vector of the case study plans
write.csv(I,file = file.path(paste0(getwd(),'/MSP_Model/Output/Data/Case_Study_Iteration_Numbers.csv')))
plans <- read.csv(file='~/Desktop/Code/MSP Planning Results April 2016/Static_plans.csv',header=F)[,I]
# Set global variables
n.sector <- 7 # Number of sectors
epsilon <- 0.2 # Stepsize of sector weights
t <- 10 # Time Horizon
r <- 0.05 # Discount rate

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
Aqua.Full.Domain.Logical <- apply(1 * (data.frame(M$Annuity,F$Annuity,K$Annuity) > 0),FUN=sum,MARGIN=1) > 0

FID <- read.csv(paste0(wkdir,'/MSP_Model/Input/Data/SeaGrant_data_complete_2015.csv'))$TARGET_FID[apply(1 * (data.frame(M$Annuity,F$Annuity,K$Annuity) > 0),FUN=sum,MARGIN=1) > 0]

# Mussel Frequency
M_plans <- plans;M_plans[M_plans != 1] <- 0
M_freq <- apply(M_plans,MARGIN = 1, FUN = sum) / length(I)

# Finfish Frequency
F_plans <- plans;F_plans[F_plans != 2] <- 0; F_plans[F_plans == 2] <- 1
F_freq <- apply(F_plans,MARGIN = 1, FUN = sum) / length(I)

# Kelp Frequency
K_plans <- plans;K_plans[K_plans != 3] <- 0; K_plans[K_plans == 3] <- 1
K_freq <- apply(K_plans,MARGIN = 1, FUN = sum) / length(I)

# Combine Data
Case_Study_Frequency <- data.frame(FID = FID, Mussel_Freq = M_freq, Finfish_Freq = F_freq, Kelp_Freq = K_freq) %>% glimpse()
write.csv(Case_Study_Frequency, file = file.path('~/Desktop/Case_Study_Frequency.csv'))
