figure2 <- function(formatList){
    source('~/MSP_Model/Scripts/Tradeoff Cartoon.r')
    figure2B <- function(format){
        if(format == 'pdf'){
            png(paste0(outfigdir,'Main_MS/pdf/Fig 2B.pdf'),width=8, height=8,units=units,res = res)
        }
        if(format == 'png'){
            png(paste0(outfigdir,'Main_MS/png/Fig 2B.png'),width=8, height=8,units=units,res = res)
        }
        if(format == 'eps'){
            postscript(paste0(outfigdir,'Main_MS/eps/Fig 2B.eps'),width=8, height=8,paper = "special")
        }
        source('~/MSP_Model/Scripts/Tradeoff Cartoon.r')
        dev.off()
    }
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
      }
      lines(sort(x.EF),y.EF[order(x.EF)],col="midnightblue",lwd=6,lty=1)
      lines(sort(x.U),y.U[order(x.U)],col = "mediumorchid1",lwd=2,lty=1)
      lines(sort(x.S),y.S[order(x.S)],col = "coral1",lwd=2,lty=1)
    }
    source(file.path(paste0(scrpdir,'pairs2.R')))
    for(format in 1:length(formatList)){
        #   pdf.options(width = 8, height = 6.4)

        color.vector=color.vector.max
        figure2B(formatList[format])
         # Color Vector For Seperating the MSP from Conventional Solutions
        # sample <- rbind(MM_test.df %>% filter(Set == 'MSP') %>% sample_n(size = 1000),
        #   MM_test.df %>% filter(Set == 'U') %>% sample_n(size = 500),
        #   MM_test.df %>% filter(Set == 'C') %>% sample_n(size = 500))
        print(formatList[format])
        if(formatList[format] == 'pdf'){
            png(paste0(outfigdir,'Main_MS/pdf/Fig 2.pdf'),width=8, height=6.4,units=units,res = res)
        }
        if(formatList[format] == 'png'){
            png(paste0(outfigdir,'Main_MS/png/Fig 2.png'),width=8, height=6.4,units=units,res = res)
        }
        if(formatList[format] == 'eps'){
            postscript(paste0(outfigdir,'Main_MS/eps/Fig 2.eps'),width=8, height=6.4,paper = "special")
        }
        pairs2(100*Master.matrix.max,lower.panel=panel.EF,
               upper.panel=NULL,col=color.vector,cex=0.8,xlim=c(0,100),
               ylim=c(0,100),pch=16,font.labels=3,cex.axis=1,las=1,xaxp=c(0,100,4),yaxp=c(0,100,2),
               gap=1)
        # title(xlab='% of Maximum',line = 1)
        title(ylab='% of Maximum')
        par(xpd=T)
        l1<-legend(.33,1,
                   legend=c('7D Frontier','2D Frontier'),fill=c("dodgerblue","midnightblue"),
                   cex=.75,title=expression(bold('Marine Spatial Planning (MSP)')),
                   title.adj = 0, bty = 'n', adj = 0, text.width=.25)
        l2<-legend(x = l1$rect$left+.0020, y = with(l1$rect, top - h)-.005,
                   legend=c('Constrained','Unconstrained'),fill=c("coral1","mediumorchid1"),
                   cex=.75,title=expression(bold('Conventional Planning ')),
                   title.adj = 0, bty = 'n', adj = 0, text.width=.25)
        inset.figure.proportion = 1/3
        inset.figure.dims = c(rep(width*(inset.figure.proportion),ts = 2))
        # WARNING!: The subplot command has changed quite drastically in recent updates to R,
        # as a result the cartoon needs to be inserted manually.
        # The next line is the now obsolete way of programmatically inserting the subplot into the tradeoff matrix.
        # try(subplot(Tradeoff.cartoon(),par = list(cex.main=2.5, cex = .45, lwd = 1)))
        par(oma=c(0,2,2,0))
        title('a', adj = 0, outer = T, cex = .75)
        title(xlab='% of Maximum',line = 3.5)
        dev.off()
        graphics.off()
    }
}
