Tradeoff <- function(sector.x = sector.x, sector.y = sector.y){
  png('~/MSP_Model/Output/Figures/Tradeoff.png')
  itor=0
  epsilon=.001
  bg = NA
  pch = 20
  cex = .08
  tolerence=10
  x=100*Master.matrix.max[[sector.x]]
  y=100*Master.matrix.max[[sector.y]]
  color.vector=color.vector.max
  x.MSP=x[color.vector=='coral']
  y.MSP=y[color.vector=='coral']
  # points(x.MSP,y.MSP, pch = 16, col = alpha("lightblue1",1/100),cex = cex)
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
  x.EF.original=x.EF;y.EF.original=y.EF;
  if(length(unique(x.EF))!=1&length(unique(y.EF))!=1){
    EF.inter=approx(x.EF.original,y=y.EF.original,n=length(alpha.mat.tmp))
    #       plot(100*EF.inter$x,100*EF.inter$y,xlab=x.name,ylab=y.name.itor,xlim=c(0,100),ylim=c(0,100),pch=19)
    x.EF=EF.inter$x;y.EF=EF.inter$y;
    #       lines(100*sort(x.EF),100*y.EF[order(x.EF)],col="black",lwd=3.5,lty=1)
    #       lines(100*sort(x.U),100*y.U[order(x.U)],col = "purple",lwd=3.5,lty=1)
    #       lines(100*sort(x.S),100*y.S[order(x.S)],col = "green",lwd=3.5,lty=3)
    #       points(100*sort(x.S),100*y.S[order(x.S)],col = "red",cex=1,pch=23,bg="black")
    #       points(100*x.EF.original,100*y.EF.original,col='red',cex=3,pch=19)
    tradeoff_exists_YN=TRUE
  }else{
    tradeoff_exists_YN=FALSE
    tradeoff.present.Y1.N0=0
    #       plot(100*x.EF.original,100*y.EF.original,xlab=x.name,ylab=y.name.itor,xlim=c(0,100),ylim=c(0,100))
    #       lines(100*sort(x.U),100*y.U[order(x.U)],col = "purple",lwd=3.5,lty=1)
    #       lines(100*sort(x.S),100*y.S[order(x.S)],col = "green",lwd=3.5,lty=3)
  }
  lwd = 1.5
  length = .25
  arr.type = 'simple'
  arr.adj = 1
  arr.length = .2
  # par(mar(c(1,4/5,4/5,2/5)))
  plot(sort(x.EF),y.EF[order(x.EF)],col="midnightblue",lwd=3.5,lty=1,xlab='% of Maximum [Finfish]',ylab='% of Maximum [Viewshed]',
       xlim=c(0,100),ylim=c(0,100),xaxp=c(0,100,4),yaxp=c(0,100,4),cex.axis=1.5,cex.lab=1.5)
  title('B', adj = 1, cex = 10)
  # points(x.MSP,y.MSP, pch = 16, col = alpha("lightblue1",1/100),cex = cex)
  lines(sort(x.U),y.U[order(x.U)],col = "mediumorchid1",lwd=3.5,lty=1)
  lines(sort(x.S),y.S[order(x.S)],col = "coral1",lwd=3.5,lty=1)

  # Unconstrained
  d.Ux=NULL
  d.Uy=NULL
  d.Ux.points=NULL
  d.Uy.points=NULL
  I.U.store=NULL
  for(itorU in 1:length(alpha.mat.tmp)){
    # Load the variables for each solution on the EF
    x.EF.U.tmp=x.EF[itorU]
    y.EF.U.tmp=y.EF[itorU]
    # Find the x.U value with the closest location to the x.EF value
    a=x.EF.U.tmp
    b=x.U
    I.U=which.min.diff(a,b)
    I.U.store[itorU]=I.U
    if(tradeoff_exists_YN==TRUE){
      d.Ux.points[itorU]=NA
      d.Uy.points[itorU]=NA
    }else{
      d.Ux.points[itorU]=x.EF.U.tmp
      d.Uy.points[itorU]=y.EF.U.tmp-unique(y.U[I.U])
    }
    I.U.store[itorU]=I.U
    if(abs(a-b[I.U])>tolerence){
      Value.of.MSP.U=NA
    }else{
      Value.of.MSP.U=y.EF.U.tmp-unique(y.U[I.U])
    }
    d.Ux[itorU]=x.EF.U.tmp
    d.Uy[itorU]=Value.of.MSP.U
  }

  d.Sx=NULL
  d.Sy=NULL
  d.Sx.points=NULL
  d.Sy.points=NULL
  I.S.store=NULL
  for(itorS in 1:length(alpha.mat.tmp)){
    # Load the variables for each solution on the EF
    x.EF.S.tmp=x.EF[itorS]
    y.EF.S.tmp=y.EF[itorS]
    # Find the x.S value with the closest location to the x.EF value
    a=x.EF.S.tmp
    b=x.S
    I.S=which.min.diff(a,b)
    if(tradeoff_exists_YN==TRUE){
      d.Sx.points[itorS]=NA
      d.Sy.points[itorS]=NA
    }else{
      d.Sy.points[itorS]=y.EF.S.tmp-unique(y.S[I.S])
      d.Sx.points[itorS]=x.EF.S.tmp
    }
    if(abs(a-b[I.S])>tolerence){
      Value.of.MSP.S=NA
    }else{
      Value.of.MSP.S=y.EF.S.tmp-unique(y.S[I.S])
      I.S.store[itorS]=I.S
    }
    d.Sx[itorS]=x.EF.S.tmp
    d.Sy[itorS]=Value.of.MSP.S
  }
  counter=0
  v0MSP.UI=which(is.na(d.Uy)==FALSE)
  v0MSP.SI=which(is.na(d.Sy)==FALSE)
  num.lines.U=seq(from=1,to=max(I.U.store),by=25)
  num.lines.S=seq(from=1,to=max(I.S.store),by=100)
  for(itor in 1:length(num.lines.U)){
    U.Value.x=x.U[num.lines.U[itor]]
    U.Value.y=y.U[num.lines.U[itor]]
  #   print(paste0('Coordinate Value Y is ',U.Value.y,' For Iteration = ',itor))
    a=U.Value.x
    b=x.EF
    I.tmp=which.min.diff(a,b)

    MSP.Value.x=x.EF[I.tmp]
    MSP.Value.y=y.EF[I.tmp]
    if(itor!=42){
      if(itor==38){
        Arrows(x0=MSP.Value.x,y0=88.712,x1=MSP.Value.x,y1=MSP.Value.y, col="grey70",lty=1,lwd=lwd, arr.type = arr.type, arr.adj = arr.adj, arr.length = arr.length)
      }else{
        Arrows(x0=MSP.Value.x,y0=U.Value.y,x1=MSP.Value.x,y1=MSP.Value.y, col="grey70",lty=1,lwd=lwd, arr.type = arr.type, arr.adj = arr.adj, arr.length = arr.length)
      }
    }

  }
  for(itor in 1:length(num.lines.S)){
    S.Value.x=x.S[num.lines.S[itor]]
    S.Value.y=y.S[num.lines.S[itor]]

    a=S.Value.x
    b=x.EF
    I.tmp=which.min.diff(a,b)

    MSP.Value.x=x.EF[I.tmp]
    MSP.Value.y=y.EF[I.tmp]
    if(MSP.Value.x>80){
      Arrows(x0=MSP.Value.x,y0=S.Value.y,x1=MSP.Value.x,y1=MSP.Value.y, col="grey60",lty=1,lwd=lwd, arr.type = arr.type, arr.adj = arr.adj, arr.length = arr.length)
    }
  }
  dev.off()
  output <- list(EF = data.frame(X = x.EF, Y = y.EF, stringsAsFactors = FALSE), U = data.frame(X = x.U, Y = y.U, stringsAsFactors = FALSE), S = data.frame(X = x.S, Y = y.S, stringsAsFactors = FALSE),
              d_U = data.frame(X = d.Ux, Y = d.Uy, stringsAsFactors = FALSE), d_S = data.frame(X = d.Sx, Y = d.Sy, stringsAsFactors = FALSE), d.U.points = data.frame(X = d.Ux.points, Y = d.Uy.points, stringsAsFactors = FALSE),
              d.S.points = data.frame(X = d.Sx.points, Y = d.Sy.points, stringsAsFactors = FALSE))
  return(output)
}
