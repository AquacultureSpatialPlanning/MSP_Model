value.of.MSP.fx_2_interpol<-function (x,y,x.name,y.name.itor,cycle_itor,epsilon=.001,tolerence=.01,...)
{
  x.MSP=x[color.vector=='coral']
  y.MSP=y[color.vector=='coral']
  x.U=x[color.vector=='purple']
  y.U=y[color.vector=='purple']
  x.S=x[color.vector=='green']
  y.S=y[color.vector=='green']
  x.EF=NULL
  y.EF=NULL
  alpha.mat.tmp=seq(from=0,by=epsilon,to=1)
  # MSP
    # Manual Version
        for(itor in 1:length(alpha.mat.tmp)){
          alpha.tmp=alpha.mat.tmp[itor]
          A=(alpha.tmp*x.MSP)+((1-alpha.tmp)*y.MSP)
          I=which(A==max(A))
          x.EF[itor]=max(unique(x.MSP[I]))
          I.tmp.x=which(x.MSP==max(unique(x.MSP[I])))
          I.tmp.y=which(y.MSP[I.tmp.x]==max(unique(y.MSP[I.tmp.x])))
          y.EF[itor]=unique(y.MSP[I.tmp.x[I.tmp.y]])}
#     # Find corrisponding alpha indices
#       unused.sector.store=names(aMatrix[,names(aMatrix)!=x.name&names(aMatrix)!=y.name.itor])
#       tmp=seq(from=1,to=nrow(aMatrix),by=1)
#       EF.I.tmp=which((aMatrix[,names(aMatrix)==x.name]>=0&aMatrix[,names(aMatrix)==y.name.itor]>=0&
#                         aMatrix[,names(aMatrix)==unused.sector.store[1]]==0&
#                         aMatrix[,names(aMatrix)==unused.sector.store[2]]==0&
#                         aMatrix[,names(aMatrix)==unused.sector.store[3]]==0&
#                         aMatrix[,names(aMatrix)==unused.sector.store[4]]==0&
#                         aMatrix[,names(aMatrix)==unused.sector.store[5]]==0))
#       EF.I=EF.I.tmp[-1]
#       x.EF=x.MSP[EF.I]
#       y.EF=y.MSP[EF.I]
  # Interpolation of the EFs
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
  # Unconstrained
    d.Ux=NULL
    d.Uy=NULL
    d.Ux.points=NULL
    d.Uy.points=NULL
    for(itorU in 1:length(alpha.mat.tmp)){
      # Load the variables for each solution on the EF
        x.EF.U.tmp=x.EF[itorU]
        y.EF.U.tmp=y.EF[itorU]
      # Find the x.U value with the closest location to the x.EF value
        a=x.EF.U.tmp
        b=x.U
        I.U=which.min.diff(a,b)
#         c.diff.U.x=(a-b)^2;
#         I.U.tmp1=which(c.diff.U.x==min(c.diff.U.x))
#         c.diff.U.y=(y.EF.U.tmp-y.U[I.U.tmp1])^2
#         I.U.tmp2=which(c.diff.U.y==min(c.diff.U.y))
#         I.U=I.U.tmp1[I.U.tmp2]
      if(tradeoff_exists_YN==TRUE){
        d.Ux.points[itorU]=NA
        d.Uy.points[itorU]=NA
      }else{
        d.Ux.points[itorU]=x.EF.U.tmp
        d.Uy.points[itorU]=y.EF.U.tmp-unique(y.U[I.U])
      }
      if(abs(a-b[I.U])>tolerence){
        Value.of.MSP.U=NA
      }else{
        Value.of.MSP.U=y.EF.U.tmp-unique(y.U[I.U])
      }
      d.Ux[itorU]=x.EF.U.tmp
      d.Uy[itorU]=Value.of.MSP.U
    }
  # Constrained
    d.Sx=NULL
    d.Sy=NULL
    d.Sx.points=NULL
    d.Sy.points=NULL
    for(itorS in 1:length(alpha.mat.tmp)){
      # Load the variables for each solution on the EF
        x.EF.S.tmp=x.EF[itorS]
        y.EF.S.tmp=y.EF[itorS]
      # Find the x.S value with the closest location to the x.EF value
        a=x.EF.S.tmp
        b=x.S
        I.S=which.min.diff(a,b)
#         c.diff.S.x=(a-b)^2;
#         I.S.tmp1=which(c.diff.S.x==min(c.diff.S.x))
#         c.diff.S.y=(y.EF.S.tmp-y.S[I.S.tmp1])^2
#         I.S.tmp2=which(c.diff.S.y==min(c.diff.S.y))
#         I.S=I.S.tmp1[I.S.tmp2]
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
      }
      d.Sx[itorS]=x.EF.S.tmp
      d.Sy[itorS]=Value.of.MSP.S
    }
  # Aquaculture Value
    Aqua.value.tmp=c(d.Ux,d.Sx)
    Aqua.value.tmp.points=c(d.Ux.points,d.Sx.points)
  # Value of MSP
    Value.of.MSP.tmp=c(d.Uy,d.Sy)
    Value.of.MSP.tmp.points=c(d.Uy.points,d.Sy.points)
  # Group (the sector on the x-axis)
    Group.tmp=rep(x.name,times=length(Aqua.value.tmp))
  # Existing sector, i.e. the sectors in which the value of MSP is measured
    Existing.sector.name.tmp=rep(y.name.itor,times=length(Aqua.value.tmp))
  # Type of Sector, i.e. aquaculture or existing, used for the faceting of variables
    type.name.tmp=y.name.itor;type.name.rep.tmp=y.name.itor
    type.name.rep.tmp[type.name.tmp=="Kelp"|type.name.tmp=="Finfish"|type.name.tmp=="Mussel"]="Aquaculture"
    type.name.rep.tmp[type.name.tmp=="Halibut"|type.name.tmp=="Viewshed"|type.name.tmp=="Benthic"|type.name.tmp=="Disease"]="Impacted Sector"
    Type.of.sector.tmp=rep(type.name.rep.tmp,times=length(Aqua.value.tmp))
  # Type of conventional
    u.string=rep('Unconstrained',times=length(d.Ux))
    s.string=rep('Constrained',times=length(d.Ux))
    Type.of.conventional=c(u.string,s.string)
  # Combine the data into a single data frame
  value.MSP.lines=data.frame(Aqua.value.tmp,Value.of.MSP.tmp,Group.tmp,Existing.sector.name.tmp,Type.of.sector.tmp,Type.of.conventional)
  value.MSP.points=data.frame(Aqua.value.tmp.points,Value.of.MSP.tmp.points,Group.tmp,Existing.sector.name.tmp,Type.of.sector.tmp,Type.of.conventional)
  value.MSP.out=list(value.MSP.lines,value.MSP.points)
  return(value.MSP.out)
}
