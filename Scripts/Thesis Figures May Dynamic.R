# This is the most concise figure script, created 3/12/16 by JS 
# Makes all five primary figures for MSP Aquaculture Paper 
# Uses maps created by Becca Gentry in GIS 
current.directory.data='C:/Users/Joel/Desktop/Thesis YTD/Code/MSP Planning Results April 2016'
# current.directory.mac='/Users/joelstevens/Desktop/December 13th work/Code/MSP Planning Results April 2016'
setwd(current.directory.data)
install.packages(c('TeachingDemos','maps','mapdata','maptools','scales','ggmap','ggplot2','grid','GGally','gridExtra','jpeg','R.matlab','png','shape','DDHFm','tiff'))
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
## Figure dimensions 
width=7
height=5.5
res=2400
units='in'
## Set Directory 
cols <- c("Mussel"="Blue","Finfish"="Salmon","Kelp"="Green","Halibut"="Burlywood","Viewshed"='Cyan',"Benthic"='Orange',"Disease"='Black') 
text.size<-12
patch.size=1.5

# Static.plans.data <- read.csv(file='Static_plans.csv',header=F)
Static.values.data <- read.csv(file='Dynamic_Values_Export.csv',header=F)
Static.plans.case.study.data <- read.csv(file='Static_plans_case_study.csv',header=F)
Static.values.case.study.data <- read.csv(file='Static_values_case_study.csv',header=F)
Static.percentage.case.study.data <- read.csv('Static_percentage_case_study.csv',header=F)
colnames(Static.values.data)=c('Mussel','Finfish','Kelp','Halibut','Viewshed','Benthic','Disease')
colnames(Static.values.case.study.data)=c('Mussel','Finfish','Kelp','Halibut','Viewshed','Benthic','Disease')
data.MSP = list(aMatrix = aMatrix.data,Static.plans = Static.plans.data,Static.plans.case.study = Static.plans.case.study.data,
              Static.values.case.study = Static.values.case.study.data,Static.percentage.case.study = Static.percentage.case.study.data)

# Conventional Models 
Unconstrained_data <- read.csv(file='Unconstrained_Dynamic_Values_April.csv',header=F)
names(Unconstrained_data)=c('Mussel','Finfish','Kelp','Halibut','Viewshed','Benthic','Disease')
Constrained_data <- read.csv(file='Constrained_Dynamic_Values_April.csv',header=F)
names(Constrained_data)=c('Mussel','Finfish','Kelp','Halibut','Viewshed','Benthic','Disease')

# Unconstrained_data <- read.csv(file='Unconstrained_Static_Values_April.csv',header=F)
# names(Unconstrained_data)=c('Mussel','Finfish','Kelp','Halibut','Viewshed','Benthic','Disease')
# Constrained_data <- read.csv(file='Constrained_Static_Values_April.csv',header=F)
# names(Constrained_data)=c('Mussel','Finfish','Kelp','Halibut','Viewshed','Benthic','Disease')

# Pure Profit Conventional Models
Unconstrained.Pure.Profit.data <- read.csv(file='Pure_Profit_Unconstrained_Dynamic_Values_April.csv',header=F)
names(Unconstrained.Pure.Profit.data)=c('Mussel','Finfish','Kelp','Halibut','Viewshed','Benthic','Disease')
Constrained.Pure.Profit.data <- read.csv(file='Pure_Profit_Constrained_Dynamic_Values_April.csv',header=F)
names(Constrained.Pure.Profit.data)=c('Mussel','Finfish','Kelp','Halibut','Viewshed','Benthic','Disease')

# Conventional and MSP
Master.matrix.max=rbind(Static.values.data,Unconstrained_data,Constrained_data)
names(Master.matrix.max)<-c('Mussel','Finfish','Kelp','Halibut','Viewshed','Benthic','Disease')
Master.matrix.max.pure.profit=rbind(Static.values.data,Unconstrained.Pure.Profit.data,Constrained.Pure.Profit.data)
color.vector.max=c(rep('coral',length.out=nrow(Static.values.data)),rep('purple',times=nrow(Unconstrained_data)),rep('green',times=nrow(Constrained_data)))


#### Plot Data 
current.directory.scripts='C:/Users/Joel/Desktop/Thesis YTD/Code'
setwd(current.directory.scripts)
## Figure 1
for(itor in 1:2){
  if(itor==1){
    png('Fig 1.png',units=units,width=width, height=height, res=res)
  }else{
    pdf('Fig 1.pdf', width=width, height=height,paper='legal')
  }
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
  
  img <- readTIFF("fig1_Stevens_v3.tif")
  g <- rasterGrob(img, interpolate=TRUE)
  
  a<-qplot(1:10, 1:10, geom="blank") + annotation_custom(g, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) + ggtitle('A') +
    theme + labs
  
  img_mussel <- readPNG("MusselValueApril.png",native=T,info=T)
  
  g_mussel <- rasterGrob(img_mussel, interpolate=TRUE)
  
  b<-qplot(1:10, 1:10, geom="blank") + ggtitle('B') + annotation_custom(g_mussel, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +
    theme + labs
  
  img_finfish <- readPNG("FishValueApril.png",native=T,info=T)
  g_finfish <- rasterGrob(img_finfish, interpolate=TRUE,just='center')
  
  c<-qplot(1:10, 1:10, geom="blank") + ggtitle('C') + annotation_custom(g_finfish, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +
    theme + labs
  
  img_kelp <- readPNG("KelpValueApril.png",native=T,info=T)
  g_kelp <- rasterGrob(img_kelp, interpolate=TRUE, just='center')
  
  d<-qplot(1:10, 1:10, geom="blank") + ggtitle('D') + annotation_custom(g_kelp, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +
    theme + labs
  grid.arrange(a,b,c,d,ncol=2,nrow=2)
  #                bottom = textGrob(expression(bold('Fig 1: ')~plain('Study domain, spatial constraints and potential value for aquaculture development. (A) Spatial constraints to aquaculture development in the Southern California Bight. (B-D) Potential value (10- year NPV) in each developable cell for mussel, finfish, and kelp aquaculture sectors.')),
  #                                  x=1,just='left'))
  dev.off() 
}
## Figure 2 
# for(itor in 2){
#   if(itor==1){
#     png('Fig 2.png',width=8, height=6.4,res=res,units=units)
#   }else{
    pdf('Fig 2.pdf',width=8, height=6.4,paper='legal')
#   }
  panel.EF<-function (x, y, itor=0, epsilon=.001, bg = NA, pch = 20, cex = .01, ...) 
  {
    x.MSP=x[color.vector=='coral']
    y.MSP=y[color.vector=='coral']
    
    points(x.MSP,y.MSP, pch = 16, col = alpha("lightblue1",1/100),cex = cex)
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
  source('pairs2.R')
  color.vector=color.vector.max # Color Vector For Seperating the MSP from Conventional Solutions
  pairs2(100*Master.matrix.max,lower.panel=panel.EF,
         upper.panel=NULL,col=color.vector,cex=0.8,xlim=c(0,100),
         ylim=c(0,100),pch=16,font.labels=3,cex.axis=1,las=1,xaxp=c(0,100,4),yaxp=c(0,100,2),
         gap=1)
  par(xpd=T)
  l1<-legend(.33,1, 
             legend=c('7D Frontier','2D Frontier'),fill=c("lightblue1","midnightblue"),
             cex=.75,title=expression(bold('Marine Spatial Planning (MSP)')),
             title.adj = 0, bty = 'n', adj = 0, text.width=.25)
  l2<-legend(x = l1$rect$left+.0020, y = with(l1$rect, top - h)-.005, 
             legend=c('Constrained','Unconstrained'),fill=c("coral1","mediumorchid1"),
             cex=.75,title=expression(bold('Conventional Management')),
             title.adj = 0, bty = 'n', adj = 0, text.width=.25)
  inset.figure.proportion = 1/3
  inset.figure.dims = c(rep(width*(inset.figure.proportion),times = 2))
  subplot(source(file='Tradeoff Cartoon.R'),x='topright',size = inset.figure.dims, type='plt', par = list(cex.main=2.5, cex = .45, lwd = 1))
  par(oma=c(0,2,2,0))
  title('A', adj = 0, outer = T, cex = .75)
  dev.off()  
# }
## Figure 3 
for(itor in 1:2){
  if(itor==1){
    png('Fig 3.png',width=width, height=height,res=res,units=units)
  }else{
    pdf('Fig 3.pdf',width=width, height=height,paper='legal')
  }
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
  
  img_hot_all <- readPNG("HotSpots_ALL_Joel.png",native=T,info=T)
  g_hot_all <- rasterGrob(img_hot_all, interpolate=TRUE)
  
  a<-qplot(1:10, 1:10, geom="blank") + ggtitle('A') +
    annotation_custom(g_hot_all, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)+
    theme + labs
  
  img_hot_mussel <- readPNG("HotSpots_Mussels_Joel.png",native=T,info=T)
  g_hot_mussel <- rasterGrob(img_hot_mussel, interpolate=TRUE)
  
  b<-qplot(1:10, 1:10, geom="blank") + ggtitle('B') +
    annotation_custom(g_hot_mussel, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)+
    theme + labs
  
  img_hot_finfish <- readPNG("HotSpots_Fish_Joel.png",native=T,info=T)
  g_hot_finfish <- rasterGrob(img_hot_finfish, interpolate=TRUE,just='center')
  
  c<-qplot(1:10, 1:10, geom="blank") + ggtitle('C') +
    annotation_custom(g_hot_finfish , xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)+
    theme + labs
  
  img_hot_kelp <- readPNG("HotSpots_Kelp_Joel_August.png",native=T,info=T)
  g_hot_kelp <- rasterGrob(img_hot_kelp, interpolate=TRUE, just='center')
  
  d<-qplot(1:10, 1:10, geom="blank") + ggtitle('D') +
    annotation_custom(g_hot_kelp, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)+
    theme + labs
  grid.arrange(a,b,c,d,ncol=2,nrow=2,padding=unit(.1,'line'))
  dev.off()
}

## Figure 4
for(itor in 1:2){
  if(itor==1){
    png('Fig 4.png',width=width, height=height,res=res,units=units)
  }else{
    pdf('Fig 4.pdf',width=width, height=height,paper='legal')
  }
  Low.impact.solutions=Static.values.case.study.data
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
  names(Low.impact.solutions)=c('ID','Values','Sector')
  Low.impact.solutions$Sector=factor(Low.impact.solutions$Sector, levels=c('Mussel','Finfish','Kelp','Halibut','Viewshed','Benthic','Disease'))
  p.bar<-ggplot(data = Low.impact.solutions,aes(x=ID,y=Values,fill=Sector,color=Sector))+
    geom_bar(stat="identity")+facet_grid(.~Sector)+ggtitle('A')+
    scale_y_continuous(labels = percent,limits=c(0,1))+
    scale_fill_manual(values=cols)+
    scale_color_manual(values=cols)+
    theme(axis.ticks.x=element_blank(),axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          panel.background=element_rect(color='white',fill='white'),panel.margin = unit(.75, "lines"),
          strip.background=element_rect(fill='white'),strip.text=element_text(size=text.size*.75),
          panel.grid=element_blank(),axis.title.y=element_text(size=text.size,color="black"),
          axis.text.y=element_text(size=text.size,color="black"),legend.title=element_text(size=text.size*.75),
          legend.text=element_text(size=text.size),plot.title=element_text(hjust=0))
  
  img_case_study_percent <- readPNG("CaseStudyPercentage_March.png",native=T,info=T)
  g_case_study_percent <- rasterGrob(img_case_study_percent, interpolate=TRUE,just='center')
  
  p.map<-qplot(1:10, 1:10, geom="blank") + ggtitle('B') +
    annotation_custom(g_case_study_percent, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)+
    theme(plot.margin = unit(c(.2,.2,.2,.2), units = "lines"),axis.text.x = element_blank(),axis.text.y = element_blank(),axis.ticks = element_blank(),
          axis.title=element_blank(),panel.background=element_rect(fill="white"),panel.grid=element_blank(),
          plot.title=element_text(hjust=0))
  
  img_case_study_plan <- readPNG("CaseStudySpecies_MarchB.png",info=T)
  g_case_study_plan <- rasterGrob(img_case_study_plan, interpolate=TRUE,just='center')
  
  p.map.case.study<-qplot(1:10, 1:10, geom="blank") + ggtitle('C') +
    annotation_custom(g_case_study_plan , xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)+
    theme(plot.margin = unit(c(.2,.2,.2,.2), units = "lines"),axis.text.x = element_blank(),axis.text.y = element_blank(),axis.ticks = element_blank(),
          axis.title=element_blank(),panel.background=element_rect(fill="white"),panel.grid=element_blank(),
          plot.title=element_text(hjust=0))
  
  p.bar.case.study <- ggplot(data = Low.impact.solutions[Low.impact.solutions$ID==576,],
                             aes(x=ID,y=Values,fill=Sector,color=Sector))+ggtitle('D')+
    geom_bar(stat="identity")+facet_grid(.~Sector)+
    scale_y_continuous(labels = percent,limits=c(0,1))+
    scale_fill_manual(values=cols)+
    scale_color_manual(values=cols)+
    theme(axis.ticks.x=element_blank(),axis.title.x=element_blank(),
          axis.title.y=element_text(size=14),axis.text.x=element_blank(),
          panel.background=element_rect(color='white',fill='white'),
          strip.background=element_rect(fill='white'),strip.text=element_blank(),
          panel.grid=element_blank(),
          axis.text.y=element_blank(),axis.title=element_blank(),axis.ticks=element_blank(),legend.position='none',
          panel.border=element_blank(),plot.title=element_text(hjust=0),plot.margin = unit(c(.2,.2,.2,.2), "cm"))
  g=ggplotGrob(p.bar.case.study) 
  p.map.case.study.full<-p.map.case.study+annotation_custom(grob=g,xmin=6.5,xmax=10,ymin=7.5)
  grid.arrange(p.bar,p.map,p.map.case.study.full,ncol=2,nrow=2,layout_matrix = rbind(c(1,1),c(2,3)))
  dev.off()
}
#   pdf.options(width = 9.5, height = 7)
## Figure 5 
# pdf.options(width = 8, height = 8)
for(itor in 1:2){
  if(itor==1){
    png('Fig 5.png',width=8, height=8,res=res,units=units)
  }else{
    pdf('Fig 5.pdf',width=8, height=8,paper='legal')
  }
  source('figure_5_code.R')
  MSP.value.data=rbind(Mussel.value.tmp,Finfish.value.tmp,Kelp.value.tmp)
  MSP.value.data.points=rbind(Mussel.value.tmp.points,Finfish.value.tmp.points,Kelp.value.tmp.points)
  names(MSP.value.data)=c('Aquaculture.Value','Value.of.MSP','Group','Sector.Name','Sector.Type','Type.of.Conventional')
  names(MSP.value.data.points)=c('Aquaculture.Value','Value.of.MSP','Group','Sector.Name','Sector.Type','Type.of.Conventional')
  MSP.value.data$Sector.Name=factor(MSP.value.data$Sector.Name,levels=c('Mussel','Finfish','Kelp','Halibut','Viewshed','Benthic','Disease'))
  Value.of.MSP.grid.plot<-ggplot(data=MSP.value.data)+
    geom_point(data=subset(MSP.value.data[as.integer(MSP.value.data$Group)!=as.integer(MSP.value.data$Sector.Name),]),
               aes(x=Aquaculture.Value,y=Value.of.MSP,shape=Type.of.Conventional,color=Type.of.Conventional),size=1.25)+
    facet_grid(Sector.Name~Group)+scale_y_continuous(labels = percent,limits=c(0,1),breaks=c(0,.5,1))+
    scale_x_continuous(labels = percent,breaks=c(0,.5,1))+
    geom_text(data=subset(MSP.value.data[as.integer(MSP.value.data$Group)==as.integer(MSP.value.data$Sector.Name),]),x=.5,y=.5,size=15,label='NA')+
    geom_point(data=MSP.value.data.points,aes(x=Aquaculture.Value,y=Value.of.MSP,shape=Type.of.Conventional,color=Type.of.Conventional),size=1.2)+
    xlab("Aquaculture Value")+ylab("Value of MSP")+
    scale_colour_manual(name = "Conventional Management:",labels=c("Constrained","Unconstrained"),values=c("coral1","mediumorchid1"))+
    scale_shape_manual(name = "Conventional Management:",labels=c("Constrained","Unconstrained"),values=c(16,24))+
    scale_fill_manual(name = "Conventional Management:",labels=c("Constrained","Unconstrained"),values=c("coral1","mediumorchid1"))+
    #     geom_line(aes(x=c(0,1),y=c(0,0)),color="grey",linetype='dashed',size=1)+
    theme_bw(base_size = 15)+theme(panel.grid=element_blank(),legend.position="bottom")
  Value.of.MSP.grid.plot+theme(axis.title=element_text(size=12),panel.margin = unit(1, "lines"),
                               strip.background=element_rect(fill='white'),strip.text=element_text(size=10),axis.ticks.margin=unit(1,'lines'))
  dev.off()
}

## SI Figures 
for(itor in 1:2){
  if(itor==1){
    png('Fig S27.png',width=8, height=6.4,res=res,units=units)
  }else{
    pdf('Fig S27.pdf',width=8, height=6.4,paper='legal')
  }
  panel.EF<-function (x, y, itor=0, epsilon=.001, bg = NA, pch = 20, cex = .01, ...) 
  {
    x.MSP=x[color.vector=='coral']
    y.MSP=y[color.vector=='coral']
    
    points(x.MSP,y.MSP, pch = 16, col = alpha("lightblue1",1/100),cex = cex)
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
  source('pairs2.R')
  color.vector=color.vector.max # Color Vector For Seperating the MSP from Conventional Solutions
  pairs2(100*Master.matrix.max.pure.profit,lower.panel=panel.EF,
         upper.panel=NULL,col=color.vector,cex=0.8,xlim=c(0,100),
         ylim=c(0,100),pch=16,font.labels=3,cex.axis=1,las=1,xaxp=c(0,100,4),yaxp=c(0,100,2),
         gap=1)
  par(xpd=T)
  l1<-legend(.33,1, 
             legend=c('7D Frontier','2D Frontier'),fill=c("lightblue1","midnightblue"),
             cex=.75,title=expression(bold('Marine Spatial Planning (MSP)')),
             title.adj = 0, bty = 'n', adj = 0, text.width=.25)
  l2<-legend(x = l1$rect$left+.0020, y = with(l1$rect, top - h)-.005, 
             legend=c('Constrained','Unconstrained')`,fill=c("coral1","mediumorchid1"),
             cex=.75,title=expression(bold('Conventional Management')),
             title.adj = 0, bty = 'n', adj = 0, text.width=.25)
  dev.off()  
}

# for(itor in 1){
#   if(itor==1){
    png('Fig S28.png',width=8, height=8,res=res,units=units)
#   }else{
#     pdf('Fig S28.pdf',width=8, height=8,paper='legal')
#   }
#   source('figure_S28_code.R')
  MSP.value.data=rbind(Mussel.value.tmp,Finfish.value.tmp,Kelp.value.tmp)
  MSP.value.data.points=rbind(Mussel.value.tmp.points,Finfish.value.tmp.points,Kelp.value.tmp.points)
  names(MSP.value.data)=c('Aquaculture.Value','Value.of.MSP','Group','Sector.Name','Sector.Type','Type.of.Conventional')
  names(MSP.value.data.points)=c('Aquaculture.Value','Value.of.MSP','Group','Sector.Name','Sector.Type','Type.of.Conventional')
  MSP.value.data$Sector.Name=factor(MSP.value.data$Sector.Name,levels=c('Mussel','Finfish','Kelp','Halibut','Viewshed','Benthic','Disease'))
  Value.of.MSP.grid.plot<-ggplot(data=MSP.value.data)+
    geom_point(data=subset(MSP.value.data[as.integer(MSP.value.data$Group)!=as.integer(MSP.value.data$Sector.Name),]),
               aes(x=Aquaculture.Value,y=Value.of.MSP,shape=Type.of.Conventional,color=Type.of.Conventional),size=1.25)+
    facet_grid(Sector.Name~Group)+scale_y_continuous(labels = percent,limits=c(0,1),breaks=c(0,.5,1))+
    scale_x_continuous(labels = percent,breaks=c(0,.5,1))+
    geom_text(data=subset(MSP.value.data[as.integer(MSP.value.data$Group)==as.integer(MSP.value.data$Sector.Name),]),x=.5,y=.5,size=15,label='NA')+
    geom_point(data=MSP.value.data.points,aes(x=Aquaculture.Value,y=Value.of.MSP,shape=Type.of.Conventional,color=Type.of.Conventional),size=1.2)+
    xlab("Aquaculture Value")+ylab("Value of MSP")+
    scale_colour_manual(name = "Conventional Management:",labels=c("Constrained","Unconstrained"),values=c("coral1","mediumorchid1"))+
    scale_shape_manual(name = "Conventional Management:",labels=c("Constrained","Unconstrained"),values=c(16,24))+
    scale_fill_manual(name = "Conventional Management:",labels=c("Constrained","Unconstrained"),values=c("coral1","mediumorchid1"))+
    #     geom_line(aes(x=c(0,1),y=c(0,0)),color="grey",linetype='dashed',size=1)+
    theme_bw(base_size = 15)+theme(panel.grid=element_blank(),legend.position="bottom")
  Value.of.MSP.grid.plot+theme(axis.title=element_text(size=12),panel.margin = unit(1, "lines"),
                               strip.background=element_rect(fill='white'),strip.text=element_text(size=10),axis.ticks.margin=unit(1,'lines'))
  dev.off()
# }






