setwd('~')
source(file.path(paste0('MSP_Model/R_Libraries.r')))
choice <- ifelse(readline(prompt = 'Install? y/n ') == 'Y',TRUE,FALSE)
print(choice)
R_Libraries(choice)
setwd('~/MSP_Model/Scripts/')

seeds <- read.csv('/Users/joel.stevens/MSP_Model/Input/Data/fid_bc_set_policies.csv')
names(seeds) <- c('FID', 'S1', 'S2', 'S3', 'S4', 'S5')
#
boo = readMat('/Users/joel.stevens/MSP_Model/Input/Data/Lester_et_al_MSPsolutions_Evals_v3/EFPayoff_a_X_wrt_DM_SeedS1S2S3.mat')
values_mat = data.frame(readMat('/Users/joel.stevens/MSP_Model/Input/Data/Raw_Impacts_FID.mat')$Raw.Impacts)

tmp = readMat('/Users/joel.stevens/MSP_Model/Input/Data/Raw_Impacts_FID.mat')
mussel = t(as.matrix(tmp[[1]]))[,2][[1]]
finfish = t(as.matrix(tmp[[1]]))[,3][[1]]
kelp = t(as.matrix(tmp[[1]]))[,4][[1]]
halibut = t(as.matrix(tmp[[1]]))[,5][[1]]
viewshed_mussel_kelp = t(as.matrix(tmp[[1]]))[,6][[1]]
viewshed_finfish = t(as.matrix(tmp[[1]]))[,7][[1]]
benthic = t(as.matrix(tmp[[1]]))[,8][[1]]
disease = t(as.matrix(tmp[[1]]))[,9][[1]]

values <- data.frame(Mussel = mussel, Finfish = finfish, Kelp = kelp, Halibut = halibut, Viewshed_Mussel = viewshed_mussel_kelp, Viewshed_Finfish = viewshed_finfish, Benthic = benthic, Disease = disease)
seeds <-

seeds <- function(values, seeds){

}







#
#
# EFPayoff_a_X_wrt_DM_SeedS1S2S3 <- setNames(data.frame(t(readMat(paste0(inpdatadir,'Lester_et_al_MSPsolutions_Evals_v3/EFPayoff_a_X_wrt_DM_SeedS1S2S3.mat'))[[1]])),
#         c('Mussel','Finfish','Kelp','Halibut','Viewshed','Benthic','Disease')) %>% mutate(Seed = c('S1','S2','S3')) %>% gather(Sector,Value,-Seed)
# p.bar <- ggplot(EFPayoff_a_X_wrt_DM_SeedS1S2S3,aes(x=factor(Sector, levels = c('Mussel','Finfish','Kelp','Halibut','Viewshed','Benthic','Disease')),y=Value,fill=factor(Seed)))+
#   geom_bar(stat="identity",position="dodge" )+
#   ylab('Value [% of maximum]') +
#   scale_fill_discrete(name="Seed",labels=c('Prioritize Existing \nSectors','Prioritize Aquaculture \nSectors','Balance Existing \nSectors and Aquaculture')) +
#   ggtitle('A') + scale_y_continuous(labels = scales::percent) +
#   theme(axis.title.x=element_blank(),
#         axis.ticks.x=element_blank(),
#         panel.background=element_rect(color='white',fill='white'),
#         strip.background=element_rect(fill='white'),strip.text=element_text(size=text.size),
#         panel.grid=element_blank(),axis.title.y=element_text(size=text.size,color="black"),
#         axis.text.x=element_text(size=text.size*.65,color="black"),axis.text.y=element_text(size=text.size*.65,color="black"),
#         legend.title=element_blank(),
#         legend.text=element_text(size=text.size * .45),
#         plot.title=element_text(hjust=0),legend.position = 'bottom',
#         legend.key.size = unit(1,"lines"),
#         legend.key = element_rect(size = 3))
# S1 <- qplot(1:10, 1:10, geom="blank") + ggtitle('B') +
#       annotation_custom(rasterGrob(readPNG(paste0(inpfigdir,"Figure5_S1.png"),native=T,info=T),interpolate = TRUE), xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)+
#       theme + labs
# S2 <- qplot(1:10, 1:10, geom="blank") + ggtitle('C') +
#       annotation_custom(rasterGrob(readPNG(paste0(inpfigdir,"Figure5_S2.png"),native=T,info=T),interpolate = TRUE), xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)+
#       theme + labs
# S3 <- qplot(1:10, 1:10, geom="blank") + ggtitle('D') +
#       annotation_custom(rasterGrob(readPNG(paste0(inpfigdir,"Figure5_S3.png"),native=T,info=T),interpolate = TRUE), xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)+
#       theme + labs
# fig5 <- arrangeGrob(p.bar,S1,S2,S3,ncol=2,nrow=2,layout_matrix = rbind(c(1,2),c(3,4)),padding=unit(-1,'line'))
# ggsave(paste0(outfigdir,'Fig 5.png'),fig5,width=width*1.15, height=height*1.15,units = units)
