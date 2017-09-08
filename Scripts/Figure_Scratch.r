# Plot Data
U.C.obj_i <- read.csv('~/MSP_Model/Input/Data/U_C_obj_i.csv',stringsAsFactors = FALSE) #%>% select(-X)
C.C.obj_i <- read.csv('~/MSP_Model/Input/Data/C_C_obj_i.csv',stringsAsFactors = FALSE) #%>% select(-X)
true_sector_names <- c('Mussel','Finfish','Kelp','Halibut','Viewshed','Benthic','Disease')
Static.values.data <- setNames(data.frame(sapply(readMat('~/MSP_Model/Input/Data/Lester_et_al_MSPsolutions_Evals_v1/EFPayoff_a_X_wrt_DM.mat'),FUN = t)),true_sector_names)
Unconstrained_data <- setNames(data.frame(sapply(readMat('~/MSP_Model/Input/Data/Lester_et_al_MSPsolutions_Evals_v1/EFPayoff_a_X_wrt_DM_UC.mat'),FUN = function(x){t(ifelse(is.na(x),0,x))})) %>%
  select(EFPayoff.a.M.wrt.DM,EFPayoff.a.F.wrt.DM,EFPayoff.a.K.wrt.DM,EFPayoff.a.H.wrt.DM,EFPayoff.a.V.wrt.DM,EFPayoff.a.B.wrt.DM,EFPayoff.a.D.wrt.DM),true_sector_names)
Constrained_data <- setNames(data.frame(sapply(readMat('~/MSP_Model/Input/Data/Lester_et_al_MSPsolutions_Evals_v1/EFPayoff_a_X_wrt_DM_CC.mat'),FUN = function(x){t(ifelse(is.na(x),0,x))})) %>%
  select(EFPayoff.a.M.wrt.DM,EFPayoff.a.F.wrt.DM,EFPayoff.a.K.wrt.DM,EFPayoff.a.H.wrt.DM,EFPayoff.a.V.wrt.DM,EFPayoff.a.B.wrt.DM,EFPayoff.a.D.wrt.DM),true_sector_names)
Master.matrix.max <- rbind(Static.values.data,Unconstrained_data,Constrained_data)
color.vector.max=c(rep('coral',length.out=nrow(Static.values.data)),rep('purple',length.out=nrow(Unconstrained_data)),rep('green',length.out=nrow(Constrained_data)))

width=7
height=5.5
res=1000
units='in'
text.size<-12
patch.size=1.5
cols <- c("Mussel"="Blue","Finfish"="Salmon","Kelp"="Green","Halibut"="Burlywood","Viewshed"='Cyan',"Benthic"='Orange',"Disease"='Black')

# current.directory.scripts=fdirs$outfigdir
# setwd(current.directory.scripts)
theme = theme(plot.margin = unit(c(.2,.2,.2,.2), units = "lines"),
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
img <- readPNG(paste0(fdirs$inpfigdir,"Fig1A Capture.png"),native=T,info=T)
g <- rasterGrob(img, interpolate=TRUE)

a<-qplot(1:10, 1:10, geom="blank") + annotation_custom(g, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) + ggtitle('A') +
  theme + labs

img_mussel <- readPNG(paste0(fdirs$inpfigdir,"MusselValueApril.png"),native=T,info=T)
g_mussel <- rasterGrob(img_mussel, interpolate=TRUE)

b <- qplot(1:10, 1:10, geom="blank") + ggtitle('B') + annotation_custom(g_mussel, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +
  theme +
  labs

img_finfish <- readPNG(paste0(fdirs$inpfigdir,"FishValueApril.png"),native=T,info=T)
g_finfish <- rasterGrob(img_finfish, interpolate=TRUE,just='center')

c<-qplot(1:10, 1:10, geom="blank") + ggtitle('C') + annotation_custom(g_finfish, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +
  theme +
  labs

img_kelp <- readPNG(paste0(fdirs$inpfigdir,"KelpValueApril.png"),native=T,info=T)
g_kelp <- rasterGrob(img_kelp, interpolate=TRUE, just='center')

d<-qplot(1:10, 1:10, geom="blank") + ggtitle('D') + annotation_custom(g_kelp, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +
  theme +
  labs

fig1 <- arrangeGrob(a,b,c,d,ncol=2,nrow=2)
# ggsave(paste0(fdirs$outfigdir,'Fig 1.png'),fig1,width=width, height=height,units = units)
ggsave(paste0(fdirs$outfigdir,'Fig 1.png'),fig1)
# Figure 2
png(paste0(fdirs$outfigdir,'Fig 2.png'),width=8, height=6.4,units=units,res = res)
panel.EF<-function (x, y, itor=0, epsilon=.001, bg = NA, pch = 20, cex = .01, ...)
{
  x.MSP=x[color.vector=='coral']
  y.MSP=y[color.vector=='coral']

  points(x.MSP,y.MSP, pch = 16, col = alpha("dodgerblue",1/75),cex = cex/2)
  x.U=x[color.vector=='purple']
  y.U=y[color.vector=='purple']
  x.S=x[color.vector=='green']
  y.S=y[color.vector=='green']
  x.EF=NULL
  y.EF=NULL
  alpha.mat.tmp=seq(from=0,by=epsilon,to=1)
  # MSP
  for(itor in 1:length(alpha.mat.tmp)){
    alpha.tmp=alpha.mat.tmp[itor]
    A=(alpha.tmp*x.MSP)+((1-alpha.tmp)*y.MSP)
    I=which(A==max(A))
    x.EF[itor]=max(unique(x.MSP[I]))
    I.tmp.x=which(x.MSP==max(unique(x.MSP[I])))
    I.tmp.y=which(y.MSP[I.tmp.x]==max(unique(y.MSP[I.tmp.x])))
    y.EF[itor]=unique(y.MSP[I.tmp.x[I.tmp.y]])}
  x.EF.original=x.EF;y.EF.original=y.EF;
  if(length(unique(x.EF.original))!=1&length(unique(x.EF.original))!=1){
    EF.inter=approx(x.EF.original,y=y.EF.original,n=length(alpha.mat.tmp))
    x.EF=EF.inter$x;y.EF=EF.inter$y;
  }else{
  }
  lines(sort(x.EF),y.EF[order(x.EF)],col="midnightblue",lwd=2,lty=1)
  lines(sort(x.U),y.U[order(x.U)],col = "mediumorchid1",lwd=2,lty=1)
  lines(sort(x.S),y.S[order(x.S)],col = "coral1",lwd=2,lty=1)
}
#   pdf.options(width = 8, height = 6.4)
source(file.path(paste0(scrpdir,'pairs2.R')))
color.vector=color.vector.max
 # Color Vector For Seperating the MSP from Conventional Solutions
# sample <- rbind(MM_test.df %>% filter(Set == 'MSP') %>% sample_n(size = 1000),
#   MM_test.df %>% filter(Set == 'U') %>% sample_n(size = 500),
#   MM_test.df %>% filter(Set == 'C') %>% sample_n(size = 500))
pairs2(100*Master.matrix.max,lower.panel=panel.EF,
       upper.panel=NULL,col=color.vector,cex=0.8,xlim=c(0,100),
       ylim=c(0,100),pch=16,font.labels=3,cex.axis=1,las=1,xaxp=c(0,100,4),yaxp=c(0,100,2),
       gap=1)
# title(xlab='% of Maximum',line = 1)
title(ylab='% of Maximum')
par(xpd=T)
l1<-legend(.33,1,
           legend=c('7D Frontier','2D Frontier'),fill=c("lightblue1","midnightblue"),
           cex=.75,title=expression(bold('Marine Spatial Planning (MSP)')),
           title.adj = 0, bty = 'n', adj = 0, text.width=.25)
l2<-legend(x = l1$rect$left+.0020, y = with(l1$rect, top - h)-.005,
           legend=c('Constrained','Unconstrained'),fill=c("coral1","mediumorchid1"),
           cex=.75,title=expression(bold('Conventional Planning ')),
           title.adj = 0, bty = 'n', adj = 0, text.width=.25)
inset.figure.proportion = 1/3
inset.figure.dims = c(rep(width*(inset.figure.proportion),ts = 2))
subplot(source(file=file.path(paste0(scrpdir,'Tradeoff Cartoonv2.R'))),x='topright',size = inset.figure.dims, type='plt', par = list(cex.main=2.5, cex = .45, lwd = 1))
par(oma=c(0,2,2,0))
title('A', adj = 0, outer = T, cex = .75)
title(xlab='% of Maximum',line = 3.5)
dev.off()
# Figure 3
# png(paste0(fdirs$outfigdir,'Fig 3.png'),width=width, height=height,res=res,units=units)
img_hot_all <- readPNG(paste0(fdirs$inpfigdir,"Figure_3_ALL.png"),native=T,info=T)
g_hot_all <- rasterGrob(img_hot_all, interpolate = TRUE)

a<-qplot(1:10, 1:10, geom="blank") + ggtitle('A') +
  annotation_custom(g_hot_all, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)+
  theme + labs

img_hot_mussel <- readPNG(paste0(fdirs$inpfigdir,"Figure_3_mussel.png"),native=T,info=T)
g_hot_mussel <- rasterGrob(img_hot_mussel, interpolate=TRUE)

b<-qplot(1:10, 1:10, geom="blank") + ggtitle('B') +
  annotation_custom(g_hot_mussel, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)+
  theme + labs

img_hot_finfish <- readPNG(paste0(fdirs$inpfigdir,"Figure_3_finfish.png"),native=T,info=T)
g_hot_finfish <- rasterGrob(img_hot_finfish, interpolate=TRUE,just='center')

c<-qplot(1:10, 1:10, geom="blank") + ggtitle('C') +
  annotation_custom(g_hot_finfish , xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)+
  theme + labs

img_hot_kelp <- readPNG(paste0(fdirs$inpfigdir,"Figure_3_kelp.png"),native=T,info=T)
g_hot_kelp <- rasterGrob(img_hot_kelp, interpolate=TRUE, just='center')

d<-qplot(1:10, 1:10, geom="blank") + ggtitle('D') +
  annotation_custom(g_hot_kelp, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)+
  theme + labs
fig4 <- arrangeGrob(a,b,c,d,ncol=2,nrow=2,padding=unit(.1,'line'))
ggsave(paste0(fdirs$outfigdir,'Fig 3.png'),fig4)
# Figure 4
EFPayoff_filter <- t(readMat('~/MSP_Model/Input/Data/Lester_et_al_MSPsolutions_Evals_v1/EFPayoff_a_X_wrt_DM_filter.mat')[[1]])
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
  scale_y_continuous(labels = percent,limits=c(0,1))+
  scale_fill_manual(values=cols)+
  scale_color_manual(values=cols)+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        panel.background=element_rect(color='white',fill='white'),panel.spacing = unit(.75, "lines"),
        strip.background=element_rect(fill='white'),strip.text=element_text(size=text.size*.75),
        panel.grid=element_blank(),axis.title.y=element_text(size=text.size,color="black"),
        axis.text.y=element_text(size=text.size,color="black"),legend.title=element_text(size=text.size*.75),
        legend.text=element_text(size=text.size),plot.title=element_text(hjust=0),plot.margin = unit(c(.2,.2,.2,.2), units = "lines"))

img_case_study_percent <- readPNG(paste0(fdirs$inpfigdir,'Figure4B.png'),native=T,info=T)
g_case_study_percent <- rasterGrob(img_case_study_percent, interpolate=TRUE,just='center')

p.map<-qplot(1:10, 1:10, geom="blank") + ggtitle('B') +
  annotation_custom(g_case_study_percent, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)+
  theme(plot.margin = unit(c(.2,1,.4,.2), units = "lines"),axis.text.x = element_blank(),axis.text.y = element_blank(),axis.ticks = element_blank(),
        axis.title=element_blank(),panel.background=element_rect(fill="white"),panel.grid=element_blank(),
        plot.title=element_text(hjust=0))
fig4 <- arrangeGrob(p.bar,p.map,ncol=1,nrow=2,layout_matrix = rbind(c(1,1),c(2,2)))
ggsave(paste0(fdirs$outfigdir,'Fig 4.png'),fig4,width=width, height=height,units = units)
# Figure 5
EFPayoff_a_X_wrt_DM_SeedS1S2S3 <- setNames(data.frame(t(readMat(paste0(fdirs$inpdatadir,'Lester_et_al_MSPsolutions_Evals_v1/EFPayoff_a_X_wrt_DM_SeedS1S2S3.mat'))[[1]])),
        c('Mussel','Finfish','Kelp','Halibut','Viewshed','Benthic','Disease')) %>% mutate(Seed = c('S1','S2','S3')) %>% gather(Sector,Value,-Seed)
p.bar <- ggplot(EFPayoff_a_X_wrt_DM_SeedS1S2S3,aes(x=factor(Sector, levels = c('Mussel','Finfish','Kelp','Halibut','Viewshed','Benthic','Disease')),y=Value,fill=factor(Seed)))+
  geom_bar(stat="identity",position="dodge")+
  scale_fill_discrete(name="Seed",labels=c('Prioritize Existing Sectors','Prioritize Aquaculture Sectors','Median Prioritization')) +
  ggtitle('A') +
  theme(axis.title.x=element_blank(),
        axis.ticks.x=element_blank(),
        panel.background=element_rect(color='white',fill='white'),panel.spacing = unit(.75, "lines"),
        strip.background=element_rect(fill='white'),strip.text=element_text(size=text.size*.75),
        panel.grid=element_blank(),axis.title.y=element_text(size=text.size,color="black"),
        axis.text.y=element_text(size=text.size,color="black"),legend.title=element_blank(),
        legend.text=element_text(size=text.size * .75),plot.title=element_text(hjust=0),legend.position = c(.17, .85))
S1 <- qplot(1:10, 1:10, geom="blank") + ggtitle('B') +
      annotation_custom(rasterGrob(readPNG(paste0(fdirs$inpfigdir,"Figure5_S1.png"),native=T,info=T),interpolate = TRUE), xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)+
      theme + labs
S2 <- qplot(1:10, 1:10, geom="blank") + ggtitle('C') +
      annotation_custom(rasterGrob(readPNG(paste0(fdirs$inpfigdir,"Figure5_S2.png"),native=T,info=T),interpolate = TRUE), xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)+
      theme + labs
S3 <- qplot(1:10, 1:10, geom="blank") + ggtitle('D') +
      annotation_custom(rasterGrob(readPNG(paste0(fdirs$inpfigdir,"Figure5_S3.png"),native=T,info=T),interpolate = TRUE), xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)+
      theme + labs
fig5 <- arrangeGrob(p.bar,S1,S2,S3,ncol=3,nrow=2,layout_matrix = rbind(c(1,1,1),c(2,3,4)))
ggsave(paste0(fdirs$outfigdir,'Fig 5.png'),fig5,width=width, height=height,units = units)
# Figure 6
source(paste0(fdirs$scrpdir,'value.of.MSP.loop_interpol.R'))
source(paste0(fdirs$scrpdir,'value.of.MSP.fx_2_interpol.R'))
source(paste0(fdirs$scrpdir,'figure_5_code.R'))
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
ggsave(paste0(fdirs$outfigdir,'Fig 6.png'),fig6,,width=8, height=8,units=units)
## Figure 1
## SI Figures
# current.directory.code <- '~/Desktop/Code/SI Figures Code/'
# current.directory.figures <- '~/Desktop/Code/SI Figures/'
# Coastline.data<-read.csv(file='Coastline.csv',header=F)
theme_3 = theme(plot.margin = unit(c(.1,.1,.1,.1), units = "lines"),
              axis.text = element_blank(),
              axis.title = element_blank(),
              axis.ticks = element_blank(),
              # axis.ticks.length = unit(0, "lines"),
              # axis.ticks.margin = unit(0, "lines"),
              panel.background=element_rect(fill="white"),
              panel.grid=element_blank(),
              plot.title=element_text(hjust =.35,vjust=.1))

theme_2 = theme(plot.margin = unit(c(.1,.1,.1,.1), units = "lines"),
              axis.text = element_blank(),
              axis.title = element_blank(),
              axis.ticks = element_blank(),
              # axis.ticks.length = unit(0, "lines"),
              # axis.ticks.margin = unit(0, "lines"),
              panel.background=element_rect(fill="white"),
              panel.grid=element_blank(),
              plot.title=element_text(hjust =.30,vjust=.1))

labs = labs(x = NULL, y = NULL)
# Figure S1
img_s1 <- readPNG(paste0(fdirs$inpfigdir,'Fig S1.png'),native=T,info=T)
g_s1 <- rasterGrob(img_s1, interpolate=TRUE)

S1 <- qplot(1:10, 1:10, geom="blank") +
      annotation_custom(g_s1, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)+
      theme_3 + labs
ggsave(paste0(fdirs$outfigdir,'Fig S1.png'),S1,width=width, height=height,units = units)
# Figure S2
img_S2A <- readPNG(paste0(fdirs$inpfigdir,'S2A.png'),native=T,info=T)
g_S2A <- rasterGrob(img_S2A, interpolate=TRUE)

S2A<-qplot(1:10, 1:10, geom="blank") + ggtitle('A') +
      annotation_custom(g_S2A, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)+
      theme_3 + labs

img_S2B <- readPNG(paste0(fdirs$inpfigdir,'S2B.png'),native=T,info=T)
g_S2B <- rasterGrob(img_S2B, interpolate=TRUE)

S2B<-qplot(1:10, 1:10, geom="blank") + ggtitle('B') +
      annotation_custom(g_S2B, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)+
      theme_3 + labs

img_S2C <- readPNG(paste0(fdirs$inpfigdir,'S2C.png'),native=T,info=T)
g_S2C <- rasterGrob(img_S2C, interpolate=TRUE)

S2C<-qplot(1:10, 1:10, geom="blank") + ggtitle('C') +
  annotation_custom(g_S2C, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)+
  theme_3 + labs

S2 <- arrangeGrob(S2A,S2B,S2C,ncol=1,nrow=3,padding=unit(-.1,'line'))
ggsave(paste0(fdirs$outfigdir,'Fig S2.png'),S2,width=width, height=height,units=units)
# S3
img_S3A <- readPNG(paste0(fdirs$inpfigdir,'S3A.png'),native=T,info=T)
g_S3A <- rasterGrob(img_S3A, interpolate=TRUE)

S3A<-qplot(1:10, 1:10, geom="blank") + ggtitle('A') +
  annotation_custom(g_S3A, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)+
  theme_3 + labs

img_S3B <- readPNG(paste0(fdirs$inpfigdir,'S3B.png'),native=T,info=T)
g_S3B <- rasterGrob(img_S3B, interpolate=TRUE)

S3B<-qplot(1:10, 1:10, geom="blank") + ggtitle('B') +
  annotation_custom(g_S3B, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)+
  theme_3 + labs

img_S3C <- readPNG(paste0(fdirs$inpfigdir,'S3C.png'),native=T,info=T)
g_S3C <- rasterGrob(img_S3C, interpolate=TRUE)

S3C<-qplot(1:10, 1:10, geom="blank") + ggtitle('C') +
  annotation_custom(g_S3C, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)+
  theme_3 + labs
S3 <- arrangeGrob(S3A,S3B,S3C,ncol=1,nrow=3,padding=unit(-.1,'line'))
ggsave(paste0(fdirs$outfigdir,'Fig S3.png'),S3,width=width, height=height,units=units)
# S4
img_S4A <- readPNG(paste0(fdirs$inpfigdir,'S5A.png'),native=T,info=T)
g_S4A <- rasterGrob(img_S4A, interpolate=TRUE)

S4A<-qplot(1:10, 1:10, geom="blank") + ggtitle('A') +
  annotation_custom(g_S4A, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)+
  theme_2 + labs

img_S4B <- readPNG(paste0(fdirs$inpfigdir,'S5B.png'),native=T,info=T)
g_S4B <- rasterGrob(img_S4B, interpolate=TRUE)

S4B<-qplot(1:10, 1:10, geom="blank") + ggtitle('B') +
  annotation_custom(g_S4B, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)+
  theme_2 + labs

S4 <- arrangeGrob(S4A,S4B,ncol=1,nrow=2,padding=unit(-.1,'line'))
ggsave(paste0(fdirs$outfigdir,'Fig S4.png'),S4,width=width, height=height,units=units)
# S5
img_S5A <- readPNG(paste0(fdirs$inpfigdir,'S6A.png'),native=T,info=T)
g_S5A <- rasterGrob(img_S5A, interpolate=TRUE)

S5A<-qplot(1:10, 1:10, geom="blank") + ggtitle('A') +
  annotation_custom(g_S5A, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)+
  theme_2 + labs

img_S5B <- readPNG(paste0(fdirs$inpfigdir,'S6B.png'),native=T,info=T)
g_S5B <- rasterGrob(img_S5B, interpolate=TRUE)

S5B<-qplot(1:10, 1:10, geom="blank") + ggtitle('B') +
  annotation_custom(g_S5B, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)+
  theme_2 + labs

S5 <- arrangeGrob(S5A,S5B,ncol=1,nrow=2,padding=unit(-.1,'line'))
ggsave(paste0(fdirs$outfigdir,'Fig S5.png'),S5,width=width, height=height,units=units)
# S6
img_S6A <- readPNG(paste0(fdirs$inpfigdir,'S7A.png'),native=T,info=T)
g_S6A <- rasterGrob(img_S6A, interpolate=TRUE)

S6A<-qplot(1:10, 1:10, geom="blank") + ggtitle('A') +
  annotation_custom(g_S6A, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)+
  theme_3 + labs + theme(plot.title=element_text(hjust=0))

img_S6B <- readPNG(paste0(fdirs$inpfigdir,'S7B.png'),native=T,info=T)
g_S6B <- rasterGrob(img_S6B, interpolate=TRUE)

S6B<-qplot(1:10, 1:10, geom="blank") + ggtitle('B') +
  annotation_custom(g_S6B, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)+
  theme_3 + labs + theme(plot.title=element_text(hjust=0))

img_S6C <- readPNG(paste0(fdirs$inpfigdir,'S7C.png'),native=T,info=T)
g_S6C <- rasterGrob(img_S6C, interpolate=TRUE)

S6C<-qplot(1:10, 1:10, geom="blank") + ggtitle('C') +
  annotation_custom(g_S6C, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)+
  theme_3 + labs + theme(plot.title=element_text(hjust=0))

img_S6D <- readPNG(paste0(fdirs$inpfigdir,'S7D.png'),native=T,info=T)
g_S6D <- rasterGrob(img_S6D, interpolate=TRUE)

S6D<-qplot(1:10, 1:10, geom="blank") + ggtitle('D') +
  annotation_custom(g_S6D, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)+
  theme_3 + labs + theme(plot.title=element_text(hjust=0))

S6 <- arrangeGrob(S6A,S6B,S6C,S6D,ncol=2,nrow=2,padding=unit(.5,'line'))
ggsave(paste0(fdirs$outfigdir,'Fig S6.png'),S6,width=width, height=height,units=units)

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
              legend.text=element_text(size=text.size),plot.title=element_text(hjust=0),legend.position = c(.17, .85))
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
              legend.text=element_text(size=text.size),plot.title=element_text(hjust=0),legend.position = c(.17, .85))
S7 <- arrangeGrob(S7A, S7B, ncol = 1, nrow = 2)
ggsave(paste0(fdirs$outfigdir,'Fig S7.png'),S7,width=7, height=7,units=units)
########################## V2
U.C.obj_i <- read.csv(file.path(paste0(inpdatadir,'U_C_obj_i.csv')),stringsAsFactors = FALSE) #%>% select(-X)
C.C.obj_i <- read.csv(file.path(paste0(inpdatadir,'C_C_obj_i.csv')),stringsAsFactors = FALSE) #%>% select(-X)
true_sector_names <- c('Mussel','Finfish','Kelp','Halibut','Viewshed','Benthic','Disease')
Static.values.data <- setNames(data.frame(sapply(readMat('~/MSP_Model/Input/Data/Lester_et_al_MSPsolutions_Evals_v3/EFPayoff_a_X_wrt_DM.mat'),FUN = t)),true_sector_names)
Unconstrained_data <- setNames(data.frame(sapply(readMat('~/MSP_Model/Output/Data/EFPayoff_a_X_wrt_DM_UC.mat'),FUN = function(x){t(ifelse(is.na(x),0,x))})) %>%
  select(EFPayoff.a.M.wrt.DM,EFPayoff.a.F.wrt.DM,EFPayoff.a.K.wrt.DM,EFPayoff.a.H.wrt.DM,EFPayoff.a.V.wrt.DM,EFPayoff.a.B.wrt.DM,EFPayoff.a.D.wrt.DM),true_sector_names)
Constrained_data <- setNames(data.frame(sapply(readMat('~/MSP_Model/Output/Data/EFPayoff_a_X_wrt_DM_CC.mat'),FUN = function(x){t(ifelse(is.na(x),0,x))})) %>%
  select(EFPayoff.a.M.wrt.DM,EFPayoff.a.F.wrt.DM,EFPayoff.a.K.wrt.DM,EFPayoff.a.H.wrt.DM,EFPayoff.a.V.wrt.DM,EFPayoff.a.B.wrt.DM,EFPayoff.a.D.wrt.DM),true_sector_names)
Master.matrix.max <- rbind(Static.values.data,Unconstrained_data,Constrained_data)
color.vector.max=c(rep('coral',length.out=nrow(Static.values.data)),rep('purple',length.out=nrow(Unconstrained_data)),rep('green',length.out=nrow(Constrained_data)))

seeds <- setNames(data.frame(t(readMat('~/MSP_Model/Input/Data/Lester_et_al_MSPsolutions_Evals_v3/EFPayoff_a_X_wrt_DM_filter.mat')[[1]])),true_sector_names) %>%
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

# current.directory.scripts=outfigdir
# setwd(current.directory.scripts)
theme = theme(plot.margin = unit(c(.2,.2,.2,.2), units = "lines"),
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
img <- readPNG(paste0(inpfigdir,"Fig1A Capture.png"),native=T,info=T)
g <- rasterGrob(img, interpolate=TRUE)

a<-qplot(1:10, 1:10, geom="blank") + annotation_custom(g, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) + ggtitle('A') +
  theme + labs

img_mussel <- readPNG(paste0(inpfigdir,"MusselValueApril.png"),native=T,info=T)
g_mussel <- rasterGrob(img_mussel, interpolate=TRUE)

b <- qplot(1:10, 1:10, geom="blank") + ggtitle('B') + annotation_custom(g_mussel, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +
  theme +
  labs

img_finfish <- readPNG(paste0(inpfigdir,"FishValueApril.png"),native=T,info=T)
g_finfish <- rasterGrob(img_finfish, interpolate=TRUE,just='center')

c<-qplot(1:10, 1:10, geom="blank") + ggtitle('C') + annotation_custom(g_finfish, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +
  theme +
  labs

img_kelp <- readPNG(paste0(inpfigdir,"KelpValueApril.png"),native=T,info=T)
g_kelp <- rasterGrob(img_kelp, interpolate=TRUE, just='center')

d<-qplot(1:10, 1:10, geom="blank") + ggtitle('D') + annotation_custom(g_kelp, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +
  theme +
  labs

fig1 <- arrangeGrob(a,b,c,d,ncol=2,nrow=2)
# ggsave(paste0(outfigdir,'Fig 1.png'),fig1,width=width, height=height,units = units)
ggsave(paste0(outfigdir,'Fig 1.png'),fig1)
# Figure 2
png(paste0(fdirs$outfigdir,'Fig 2.png'),width=8, height=6.4,units=units,res = res)
panel.EF<-function (x, y, itor=0, epsilon=.001, bg = NA, pch = 20, cex = .01, ...)
{
  x.MSP=x[color.vector=='coral']
  y.MSP=y[color.vector=='coral']

  points(x.MSP,y.MSP, pch = 16, col = alpha("dodgerblue",1/75),cex = cex/2)
  x.U=x[color.vector=='purple']
  y.U=y[color.vector=='purple']
  x.S=x[color.vector=='green']
  y.S=y[color.vector=='green']
  x.EF=NULL
  y.EF=NULL
  alpha.mat.tmp=seq(from=0,by=epsilon,to=1)
  # MSP
  for(itor in 1:length(alpha.mat.tmp)){
    alpha.tmp=alpha.mat.tmp[itor]
    A=(alpha.tmp*x.MSP)+((1-alpha.tmp)*y.MSP)
    I=which(A==max(A))
    x.EF[itor]=max(unique(x.MSP[I]))
    I.tmp.x=which(x.MSP==max(unique(x.MSP[I])))
    I.tmp.y=which(y.MSP[I.tmp.x]==max(unique(y.MSP[I.tmp.x])))
    y.EF[itor]=unique(y.MSP[I.tmp.x[I.tmp.y]])}
  x.EF.original=x.EF;y.EF.original=y.EF;
  if(length(unique(x.EF.original))!=1&length(unique(x.EF.original))!=1){
    EF.inter=approx(x.EF.original,y=y.EF.original,n=length(alpha.mat.tmp))
    x.EF=EF.inter$x;y.EF=EF.inter$y;
  }else{
  }
  lines(sort(x.EF),y.EF[order(x.EF)],col="midnightblue",lwd=2,lty=1)
  lines(sort(x.U),y.U[order(x.U)],col = "mediumorchid1",lwd=2,lty=1)
  lines(sort(x.S),y.S[order(x.S)],col = "coral1",lwd=2,lty=1)
}
#   pdf.options(width = 8, height = 6.4)
source(file.path(paste0(scrpdir,'pairs2.R')))
color.vector=color.vector.max
 # Color Vector For Seperating the MSP from Conventional Solutions
# sample <- rbind(MM_test.df %>% filter(Set == 'MSP') %>% sample_n(size = 1000),
#   MM_test.df %>% filter(Set == 'U') %>% sample_n(size = 500),
#   MM_test.df %>% filter(Set == 'C') %>% sample_n(size = 500))
pairs2(100*Master.matrix.max,lower.panel=panel.EF,
       upper.panel=NULL,col=color.vector,cex=0.8,xlim=c(0,100),
       ylim=c(0,100),pch=16,font.labels=3,cex.axis=1,las=1,xaxp=c(0,100,4),yaxp=c(0,100,2),
       gap=1)
# title(xlab='% of Maximum',line = 1)
title(ylab='% of Maximum')
par(xpd=T)
l1<-legend(.33,1,
           legend=c('7D Frontier','2D Frontier'),fill=c("lightblue1","midnightblue"),
           cex=.75,title=expression(bold('Marine Spatial Planning (MSP)')),
           title.adj = 0, bty = 'n', adj = 0, text.width=.25)
l2<-legend(x = l1$rect$left+.0020, y = with(l1$rect, top - h)-.005,
           legend=c('Constrained','Unconstrained'),fill=c("coral1","mediumorchid1"),
           cex=.75,title=expression(bold('Conventional Planning ')),
           title.adj = 0, bty = 'n', adj = 0, text.width=.25)
inset.figure.proportion = 1/3
inset.figure.dims = c(rep(width*(inset.figure.proportion),ts = 2))
subplot(source(file=file.path(paste0(scrpdir,'Tradeoff Cartoonv2.R'))),x='topright',size = inset.figure.dims, type='plt', par = list(cex.main=2.5, cex = .45, lwd = 1))
par(oma=c(0,2,2,0))
title('A', adj = 0, outer = T, cex = .75)
title(xlab='% of Maximum',line = 3.5)
dev.off()
# Figure 3
# Hot spot set up
system2(matlab_root,
  args = c('-nodesktop','-noFigureWindows','-nosplash','-r',
  paste0("run\\(\\'",scrpdir,"Scratch_File.m\\'\\)")))
source(file.path(paste0(scrpdir,'Scratch_File.r')))
# png(paste0(outfigdir,'Fig 3.png'),width=width, height=height,res=res,units=units)
img_hot_all <- readPNG(paste0(inpfigdir,"Figure_3_ALL.png"),native=T,info=T)
g_hot_all <- rasterGrob(img_hot_all, interpolate = TRUE)

a<-qplot(1:10, 1:10, geom="blank") + ggtitle('A') +
  annotation_custom(g_hot_all, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)+
  theme + labs

img_hot_mussel <- readPNG(paste0(inpfigdir,"Figure_3_mussel.png"),native=T,info=T)
g_hot_mussel <- rasterGrob(img_hot_mussel, interpolate=TRUE)

b<-qplot(1:10, 1:10, geom="blank") + ggtitle('B') +
  annotation_custom(g_hot_mussel, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)+
  theme + labs

img_hot_finfish <- readPNG(paste0(inpfigdir,"Figure_3_finfish.png"),native=T,info=T)
g_hot_finfish <- rasterGrob(img_hot_finfish, interpolate=TRUE,just='center')

c<-qplot(1:10, 1:10, geom="blank") + ggtitle('C') +
  annotation_custom(g_hot_finfish , xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)+
  theme + labs

img_hot_kelp <- readPNG(paste0(inpfigdir,"Figure_3_kelp.png"),native=T,info=T)
g_hot_kelp <- rasterGrob(img_hot_kelp, interpolate=TRUE, just='center')

d<-qplot(1:10, 1:10, geom="blank") + ggtitle('D') +
  annotation_custom(g_hot_kelp, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)+
  theme + labs
fig4 <- arrangeGrob(a,b,c,d,ncol=2,nrow=2,padding=unit(.1,'line'))
ggsave(paste0(outfigdir,'Fig 3.png'),fig4)
# Figure 4
EFPayoff_filter <- t(readMat('~/MSP_Model/Input/Data/Lester_et_al_MSPsolutions_Evals_v3/EFPayoff_a_X_wrt_DM_filter.mat')[[1]])
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
  scale_y_continuous(labels = percent,limits=c(0,1))+
  scale_fill_manual(values=cols)+
  scale_color_manual(values=cols)+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        panel.background=element_rect(color='white',fill='white'),panel.spacing = unit(0, "lines"),
        strip.background=element_rect(fill='white'),strip.text=element_text(size=text.size*.75),
        panel.grid=element_blank(),axis.title.y=element_text(size=text.size,color="black"),
        axis.text.y=element_text(size=text.size,color="black"),legend.title=element_text(size=text.size),
        legend.text=element_text(size=text.size),plot.title=element_text(hjust=0),plot.margin = unit(c(.2,.2,.2,.2), units = "lines"))

img_case_study_percent <- readPNG(paste0(inpfigdir,'Figure4B.png'),native=T,info=T)
g_case_study_percent <- rasterGrob(img_case_study_percent, interpolate=TRUE,just='center')

p.map<-qplot(1:10, 1:10, geom="blank") + ggtitle('B') +
  annotation_custom(g_case_study_percent, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)+
  theme(plot.margin = unit(c(.2,1,.4,.2), units = "lines"),axis.text.x = element_blank(),axis.text.y = element_blank(),axis.ticks = element_blank(),
        axis.title=element_blank(),panel.background=element_rect(fill="white"),panel.grid=element_blank(),
        plot.title=element_text(hjust=0))
fig4 <- arrangeGrob(p.bar,p.map,ncol=1,nrow=2,layout_matrix = rbind(c(1,1),c(2,2)))
ggsave(paste0(outfigdir,'Fig 4.png'),fig4,width=width, height=height,units = units)
# Figure 5
EFPayoff_a_X_wrt_DM_SeedS1S2S3 <- setNames(data.frame(t(readMat(paste0(inpdatadir,'Lester_et_al_MSPsolutions_Evals_v3/EFPayoff_a_X_wrt_DM_SeedS1S2S3.mat'))[[1]])),
        c('Mussel','Finfish','Kelp','Halibut','Viewshed','Benthic','Disease')) %>% mutate(Seed = c('S1','S2','S3')) %>% gather(Sector,Value,-Seed)
p.bar <- ggplot(EFPayoff_a_X_wrt_DM_SeedS1S2S3,aes(x=factor(Sector, levels = c('Mussel','Finfish','Kelp','Halibut','Viewshed','Benthic','Disease')),y=Value,fill=factor(Seed)))+
  geom_bar(stat="identity",position="dodge")+
  scale_fill_discrete(name="Seed",labels=c('Prioritize Existing Sectors','Prioritize Aquaculture Sectors','Balance Existing Sectors and Aquaculture')) +
  ggtitle('A') +
  theme(axis.title.x=element_blank(),
        axis.ticks.x=element_blank(),
        panel.background=element_rect(color='white',fill='white'),panel.spacing = unit(.75, "lines"),
        strip.background=element_rect(fill='white'),strip.text=element_text(size=text.size*.75),
        panel.grid=element_blank(),axis.title.y=element_text(size=text.size,color="black"),
        axis.text.y=element_text(size=text.size,color="black"),legend.title=element_blank(),
        legend.text=element_text(size=text.size * .75),plot.title=element_text(hjust=0),legend.position = c(.17, .85))
S1 <- qplot(1:10, 1:10, geom="blank") + ggtitle('B') +
      annotation_custom(rasterGrob(readPNG(paste0(inpfigdir,"Figure5_S1.png"),native=T,info=T),interpolate = TRUE), xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)+
      theme + labs
S2 <- qplot(1:10, 1:10, geom="blank") + ggtitle('C') +
      annotation_custom(rasterGrob(readPNG(paste0(inpfigdir,"Figure5_S2.png"),native=T,info=T),interpolate = TRUE), xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)+
      theme + labs
S3 <- qplot(1:10, 1:10, geom="blank") + ggtitle('D') +
      annotation_custom(rasterGrob(readPNG(paste0(inpfigdir,"Figure5_S3.png"),native=T,info=T),interpolate = TRUE), xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)+
      theme + labs
fig5 <- arrangeGrob(p.bar,S1,S2,S3,ncol=2,nrow=2,layout_matrix = rbind(c(1,2),c(3,4)))
ggsave(paste0(outfigdir,'Fig 5.png'),fig5,width=width, height=height,units = units)
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
ggsave(paste0(outfigdir,'Fig 6.png'),fig6,width=8, height=8,units=units)
## Figure 1
## SI Figures
# current.directory.code <- '~/Desktop/Code/SI Figures Code/'
# current.directory.figures <- '~/Desktop/Code/SI Figures/'
# Coastline.data<-read.csv(file='Coastline.csv',header=F)
theme_3 = theme(plot.margin = unit(c(.1,.1,.1,.1), units = "lines"),
              axis.text = element_blank(),
              axis.title = element_blank(),
              axis.ticks = element_blank(),
              # axis.ticks.length = unit(0, "lines"),
              # axis.ticks.margin = unit(0, "lines"),
              panel.background=element_rect(fill="white"),
              panel.grid=element_blank(),
              plot.title=element_text(hjust =.35,vjust=.1))

theme_2 = theme(plot.margin = unit(c(.1,.1,.1,.1), units = "lines"),
              axis.text = element_blank(),
              axis.title = element_blank(),
              axis.ticks = element_blank(),
              # axis.ticks.length = unit(0, "lines"),
              # axis.ticks.margin = unit(0, "lines"),
              panel.background=element_rect(fill="white"),
              panel.grid=element_blank(),
              plot.title=element_text(hjust =.30,vjust=.1))

labs = labs(x = NULL, y = NULL)
# Figure S1
img_s1 <- readPNG(paste0(inpfigdir,'Fig S1.png'),native=T,info=T)
g_s1 <- rasterGrob(img_s1, interpolate=TRUE)

S1 <- qplot(1:10, 1:10, geom="blank") +
      annotation_custom(g_s1, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)+
      theme_3 + labs
ggsave(paste0(outfigdir,'Fig S1.png'),S1,width=width, height=height,units = units)
# Figure S2
img_S2A <- readPNG(paste0(inpfigdir,'S2A.png'),native=T,info=T)
g_S2A <- rasterGrob(img_S2A, interpolate=TRUE)

S2A<-qplot(1:10, 1:10, geom="blank") + ggtitle('A') +
      annotation_custom(g_S2A, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)+
      theme_3 + labs

img_S2B <- readPNG(paste0(inpfigdir,'S2B.png'),native=T,info=T)
g_S2B <- rasterGrob(img_S2B, interpolate=TRUE)

S2B<-qplot(1:10, 1:10, geom="blank") + ggtitle('B') +
      annotation_custom(g_S2B, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)+
      theme_3 + labs

img_S2C <- readPNG(paste0(inpfigdir,'S2C.png'),native=T,info=T)
g_S2C <- rasterGrob(img_S2C, interpolate=TRUE)

S2C<-qplot(1:10, 1:10, geom="blank") + ggtitle('C') +
  annotation_custom(g_S2C, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)+
  theme_3 + labs

S2 <- arrangeGrob(S2A,S2B,S2C,ncol=1,nrow=3,padding=unit(-.1,'line'))
ggsave(paste0(outfigdir,'Fig S2.png'),S2,width=width, height=height,units=units)
# S3
img_S3A <- readPNG(paste0(inpfigdir,'S3A.png'),native=T,info=T)
g_S3A <- rasterGrob(img_S3A, interpolate=TRUE)

S3A<-qplot(1:10, 1:10, geom="blank") + ggtitle('A') +
  annotation_custom(g_S3A, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)+
  theme_3 + labs

img_S3B <- readPNG(paste0(inpfigdir,'S3B.png'),native=T,info=T)
g_S3B <- rasterGrob(img_S3B, interpolate=TRUE)

S3B<-qplot(1:10, 1:10, geom="blank") + ggtitle('B') +
  annotation_custom(g_S3B, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)+
  theme_3 + labs

img_S3C <- readPNG(paste0(inpfigdir,'S3C.png'),native=T,info=T)
g_S3C <- rasterGrob(img_S3C, interpolate=TRUE)

S3C<-qplot(1:10, 1:10, geom="blank") + ggtitle('C') +
  annotation_custom(g_S3C, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)+
  theme_3 + labs
S3 <- arrangeGrob(S3A,S3B,S3C,ncol=1,nrow=3,padding=unit(-.1,'line'))
ggsave(paste0(outfigdir,'Fig S3.png'),S3,width=width, height=height,units=units)
# S4
img_S4A <- readPNG(paste0(inpfigdir,'S5A.png'),native=T,info=T)
g_S4A <- rasterGrob(img_S4A, interpolate=TRUE)

S4A<-qplot(1:10, 1:10, geom="blank") + ggtitle('A') +
  annotation_custom(g_S4A, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)+
  theme_2 + labs

img_S4B <- readPNG(paste0(inpfigdir,'S5B.png'),native=T,info=T)
g_S4B <- rasterGrob(img_S4B, interpolate=TRUE)

S4B<-qplot(1:10, 1:10, geom="blank") + ggtitle('B') +
  annotation_custom(g_S4B, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)+
  theme_2 + labs

S4 <- arrangeGrob(S4A,S4B,ncol=1,nrow=2,padding=unit(-.1,'line'))
ggsave(paste0(outfigdir,'Fig S4.png'),S4,width=width, height=height,units=units)
# S5
img_S5A <- readPNG(paste0(inpfigdir,'S6A.png'),native=T,info=T)
g_S5A <- rasterGrob(img_S5A, interpolate=TRUE)

S5A<-qplot(1:10, 1:10, geom="blank") + ggtitle('A') +
  annotation_custom(g_S5A, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)+
  theme_2 + labs

img_S5B <- readPNG(paste0(inpfigdir,'S6B.png'),native=T,info=T)
g_S5B <- rasterGrob(img_S5B, interpolate=TRUE)

S5B<-qplot(1:10, 1:10, geom="blank") + ggtitle('B') +
  annotation_custom(g_S5B, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)+
  theme_2 + labs

S5 <- arrangeGrob(S5A,S5B,ncol=1,nrow=2,padding=unit(-.1,'line'))
ggsave(paste0(outfigdir,'Fig S5.png'),S5,width=width, height=height,units=units)
# S6
img_S6A <- readPNG(paste0(inpfigdir,'S7A.png'),native=T,info=T)
g_S6A <- rasterGrob(img_S6A, interpolate=TRUE)

S6A<-qplot(1:10, 1:10, geom="blank") + ggtitle('A') +
  annotation_custom(g_S6A, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)+
  theme_3 + labs + theme(plot.title=element_text(hjust=0))

img_S6B <- readPNG(paste0(inpfigdir,'S7B.png'),native=T,info=T)
g_S6B <- rasterGrob(img_S6B, interpolate=TRUE)

S6B<-qplot(1:10, 1:10, geom="blank") + ggtitle('B') +
  annotation_custom(g_S6B, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)+
  theme_3 + labs + theme(plot.title=element_text(hjust=0))

img_S6C <- readPNG(paste0(inpfigdir,'S7C.png'),native=T,info=T)
g_S6C <- rasterGrob(img_S6C, interpolate=TRUE)

S6C<-qplot(1:10, 1:10, geom="blank") + ggtitle('C') +
  annotation_custom(g_S6C, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)+
  theme_3 + labs + theme(plot.title=element_text(hjust=0))

img_S6D <- readPNG(paste0(inpfigdir,'S7D.png'),native=T,info=T)
g_S6D <- rasterGrob(img_S6D, interpolate=TRUE)

S6D<-qplot(1:10, 1:10, geom="blank") + ggtitle('D') +
  annotation_custom(g_S6D, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)+
  theme_3 + labs + theme(plot.title=element_text(hjust=0))

S6 <- arrangeGrob(S6A,S6B,S6C,S6D,ncol=2,nrow=2,padding=unit(.5,'line'))
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
              legend.text=element_text(size=text.size),plot.title=element_text(hjust=0),legend.position = c(.17, .85))
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
              legend.text=element_text(size=text.size),plot.title=element_text(hjust=0),legend.position = c(.17, .85))
S7 <- arrangeGrob(S7A, S7B, ncol = 1, nrow = 2)
ggsave(paste0(outfigdir,'Fig S7.png'),S7,width=7, height=7,units=units)
