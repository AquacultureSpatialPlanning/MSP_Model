# Constrained conventional model
C.df.tmp <- C.df %>%
            mutate(FID = fulldomain[Aqua.Full.Domain.Logical]) %>%
            select(FID, M, F, K) %>% glimpse()
C.C.Full.df <- lapply(1:3,FUN = function(x){
  setNames(data.frame(cbind(1:i,fulldomain[Aqua.Full.Domain.Logical],C.df[,x])),c('i','FID','Value')) %>% mutate(Aqua_Type = c('Mussel','Finfish','Kelp')[x]) %>% arrange(desc(Value))
  # filter(setNames(data.frame(cbind(1:i,fulldomain[Aqua.Full.Domain.Logical],C.df[,x])),c('i','FID','Value')),Value > 0) %>% arrange(desc(Value))
  })
# constrained_c <- function(x){
#   print(x)
#   num_sites <- x # Number of sites to be developed in the spatial plan
#   C.C.Full.df.tmp <- C.C.Full.df
#   tmp <- rep(0,i) # empty spatial plan
#   num_dev <- 0
#   repeat{
#     df <- filter(data.frame(rbindlist(lapply(C.C.Full.df.tmp,'[',1,))),Value > 0) # Make a data frame from the most suitable site of the three aquaculture forms, remove aqua sectors which are zero
#     # print(nrow(df))
#     current_choices <- c('Mussel','Finfish','Kelp')
#     for(itor in 1:nrow(df)){
#       df.dev <- df[which.max(df$Value),] # The first sector developed
#       df.dev.FID <- df.dev$FID # The FID of the first aquaculture site to be developed in the triplet
#       choice <- ifelse(df.dev$Aqua_Type != 'Mussel',ifelse(df.dev$Aqua_Type != 'Kelp',2,3),1)
#       current_choices <- current_choices[which(current_choices != df.dev$Aqua_Type)]
#       tmp[which(fulldomain[Aqua.Full.Domain.Logical] == df.dev.FID)] <- choice # Develop that site for the selected type of aquaculture
#       # print(table(tmp))
#       # Remove that FID from the ordered suitability for all three aqauculture types
#       C.C.Full.df.tmp <- lapply(C.C.Full.df.tmp,FUN = function(y){arrange(filter(y,FID != df.dev.FID),desc(Value))})
#       df <- filter(filter(data.frame(rbindlist(lapply(C.C.Full.df.tmp,'[',1,))),Value > 0),FID != df.dev.FID & Aqua_Type %in% current_choices) # Update the temporary dataframe for the current triplet
#       if(length(which(tmp != 0)) == num_sites){break}
#     }
#     if(length(which(tmp != 0)) == num_sites){break}
#   }
#   spatial_plan <- tmp
#   print(table(spatial_plan))
#   return(spatial_plan)
# }
# constrained_c(1061)
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
