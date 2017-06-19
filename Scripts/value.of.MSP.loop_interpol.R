value.of.MSP.loop_interpol<-function (x.name,data.tmp,...) 
{
  x.value=data.Value.of.MSP.Ratio[,names(Static.values)==x.name]
  x=x.value
  y.names=names(data.Value.of.MSP.Ratio[,names(Static.values)!=x.name])
  y.values=data.Value.of.MSP.Ratio[,names(Static.values)!=x.name]
  for(itor_cycle in 1:length(y.names)){
    y=y.values[,itor_cycle]
    y.name.itor=y.names[itor_cycle]
    # setwd("C:/Users/Joel/Desktop/Thesis Master/Models/Current Models/Ratio Conventional/R Files")#<-Change directory 
    value.MSP.out<-value.of.MSP.fx_2_interpol(x,y,x.name,y.name.itor,cycle_itor,epsilon=.001,tolerence=.01)
    nam1 <- paste("Store.lines.tmp", itor_cycle, sep = "")
    assign(nam1, value.MSP.out[[1]])
    nam2 <- paste("Store.points.tmp", itor_cycle, sep = "")
    assign(nam2, value.MSP.out[[2]])
  }
  loopOut.lines=rbind(Store.lines.tmp1,Store.lines.tmp2,Store.lines.tmp3,Store.lines.tmp4,Store.lines.tmp5,Store.lines.tmp6)
  loopOut.points=rbind(Store.points.tmp1,Store.points.tmp2,Store.points.tmp3,Store.points.tmp4,Store.points.tmp5,Store.points.tmp6)
  loopOut=list(loopOut.lines,loopOut.points)
  return(loopOut)
}
