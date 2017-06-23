Static.values=Static.values.data
color.vector=color.vector.max
data.Value.of.MSP.Ratio.sum.impacts=Master.matrix.max
data.Value.of.MSP.Ratio=Master.matrix.max
for(sector_itor in 1:3){
  if(sector_itor==1){
    #Mussel is selected to calculate the value of MSP 
    x.name="Mussel"
    data.tmp=data.Value.of.MSP.Ratio
    loopOut<-value.of.MSP.loop_interpol(x.name,data.tmp)
    Mussel.value.tmp=loopOut[[1]] 
    Mussel.value.tmp.points=loopOut[[2]] 
    # For the X's
    d.Ux=seq(from=0,to=1,length.out=nrow(Mussel.value.tmp[Mussel.value.tmp$Existing.sector.name.tmp=='Finfish'&Mussel.value.tmp$Type.of.conventional=='Constrained',]))
    d.Uy=seq(from=1,to=0,length.out=nrow(Mussel.value.tmp[Mussel.value.tmp$Existing.sector.name.tmp=='Finfish'&Mussel.value.tmp$Type.of.conventional=='Constrained',]))
    d.Sx=seq(from=0,to=1,length.out=nrow(Mussel.value.tmp[Mussel.value.tmp$Existing.sector.name.tmp=='Finfish'&Mussel.value.tmp$Type.of.conventional=='Constrained',]))
    d.Sy=seq(from=0,to=1,length.out=nrow(Mussel.value.tmp[Mussel.value.tmp$Existing.sector.name.tmp=='Finfish'&Mussel.value.tmp$Type.of.conventional=='Constrained',]))
    # Aquaculture Value 
    Aqua.value.tmp=c(d.Ux,d.Sx)
    # Value of MSP
    Value.of.MSP.tmp=c(d.Uy,d.Sy)
    # Group (the sector on the x-axis)
    Group.tmp=rep(x.name,times=length(Aqua.value.tmp))
    # Existing sector, i.e. the sectors in which the value of MSP is measured 
    Existing.sector.name.tmp=rep(x.name,times=length(Aqua.value.tmp))
    # Type of Sector, i.e. aquaculture or existing, used for the faceting of variables 
    type.name.tmp=x.name;type.name.rep.tmp=x.name
    type.name.rep.tmp[type.name.tmp=="Kelp"|type.name.tmp=="Finfish"|type.name.tmp=="Mussel"]="Aquaculture"
    type.name.rep.tmp[type.name.tmp=="Halibut"|type.name.tmp=="Viewshed"|type.name.tmp=="Benthic"|type.name.tmp=="Disease"]="Impacted Sector"
    Type.of.sector.tmp=rep(type.name.rep.tmp,times=length(Aqua.value.tmp))
    # Type of conventional 
    u.string=rep('Unconstrained',times=length(d.Ux))
    s.string=rep('Unconstrained',times=length(d.Ux))
    Type.of.conventional=c(u.string,s.string)
    value.MSP.out.mussel=data.frame(Aqua.value.tmp,Value.of.MSP.tmp,Group.tmp,Existing.sector.name.tmp,Type.of.sector.tmp,Type.of.conventional)
    Mussel.value.tmp=rbind(Mussel.value.tmp,value.MSP.out.mussel)
  }
  if(sector_itor==2){
    #Finfish is selected to calculate the value of MSP 
    x.name="Finfish"
    data.tmp=data.Value.of.MSP.Ratio
    loopOut<-value.of.MSP.loop_interpol(x.name,data.tmp)
    Finfish.value.tmp=loopOut[[1]]
    Finfish.value.tmp.points=loopOut[[2]]
    # For the X's
    d.Ux=seq(from=0,to=1,length.out=nrow(Finfish.value.tmp[Finfish.value.tmp$Existing.sector.name.tmp=='Kelp'&Finfish.value.tmp$Type.of.conventional=='Constrained',]))
    d.Uy=seq(from=1,to=0,length.out=nrow(Finfish.value.tmp[Finfish.value.tmp$Existing.sector.name.tmp=='Kelp'&Finfish.value.tmp$Type.of.conventional=='Constrained',]))
    d.Sx=seq(from=0,to=1,length.out=nrow(Finfish.value.tmp[Finfish.value.tmp$Existing.sector.name.tmp=='Kelp'&Finfish.value.tmp$Type.of.conventional=='Constrained',]))
    d.Sy=seq(from=0,to=1,length.out=nrow(Finfish.value.tmp[Finfish.value.tmp$Existing.sector.name.tmp=='Kelp'&Finfish.value.tmp$Type.of.conventional=='Constrained',]))
    # Aquaculture Value 
    Aqua.value.tmp=c(d.Ux,d.Sx)
    # Value of MSP
    Value.of.MSP.tmp=c(d.Uy,d.Sy)
    # Group (the sector on the x-axis)
    Group.tmp=rep(x.name,times=length(Aqua.value.tmp))
    # Existing sector, i.e. the sectors in which the value of MSP is measured 
    Existing.sector.name.tmp=rep(x.name,times=length(Aqua.value.tmp))
    # Type of Sector, i.e. aquaculture or existing, used for the faceting of variables 
    type.name.tmp=x.name;type.name.rep.tmp=x.name
    type.name.rep.tmp[type.name.tmp=="Kelp"|type.name.tmp=="Finfish"|type.name.tmp=="Mussel"]="Aquaculture"
    type.name.rep.tmp[type.name.tmp=="Halibut"|type.name.tmp=="Viewshed"|type.name.tmp=="Benthic"|type.name.tmp=="Disease"]="Impacted Sector"
    Type.of.sector.tmp=rep(type.name.rep.tmp,times=length(Aqua.value.tmp))
    # Type of conventional 
    u.string=rep('Unconstrained',times=length(d.Ux))
    s.string=rep('Unconstrained',times=length(d.Ux))
    Type.of.conventional=c(u.string,s.string)
    value.MSP.out.finfish=data.frame(Aqua.value.tmp,Value.of.MSP.tmp,Group.tmp,Existing.sector.name.tmp,Type.of.sector.tmp,Type.of.conventional)
    Finfish.value.tmp=rbind(Finfish.value.tmp,value.MSP.out.finfish)
  }
  if(sector_itor==3){
    #Kelp is selected to calculate the value of MSP 
    x.name="Kelp"
    data.tmp=data.Value.of.MSP.Ratio
    loopOut<-value.of.MSP.loop_interpol(x.name,data.tmp)
    Kelp.value.tmp=loopOut[[1]]
    Kelp.value.tmp.points=loopOut[[2]]
    # For the X's
    d.Ux=seq(from=0,to=1,length.out=nrow(Kelp.value.tmp[Kelp.value.tmp$Existing.sector.name.tmp=='Mussel'&Kelp.value.tmp$Type.of.conventional=='Constrained',]))
    d.Uy=seq(from=1,to=0,length.out=nrow(Kelp.value.tmp[Kelp.value.tmp$Existing.sector.name.tmp=='Mussel'&Kelp.value.tmp$Type.of.conventional=='Constrained',]))
    d.Sx=seq(from=0,to=1,length.out=nrow(Kelp.value.tmp[Kelp.value.tmp$Existing.sector.name.tmp=='Mussel'&Kelp.value.tmp$Type.of.conventional=='Constrained',]))
    d.Sy=seq(from=0,to=1,length.out=nrow(Kelp.value.tmp[Kelp.value.tmp$Existing.sector.name.tmp=='Mussel'&Kelp.value.tmp$Type.of.conventional=='Constrained',]))
    # Aquaculture Value 
    Aqua.value.tmp=c(d.Ux,d.Sx)
    # Value of MSP
    Value.of.MSP.tmp=c(d.Uy,d.Sy)
    # Group (the sector on the x-axis)
    Group.tmp=rep(x.name,times=length(Aqua.value.tmp))
    # Existing sector, i.e. the sectors in which the value of MSP is measured 
    Existing.sector.name.tmp=rep(x.name,times=length(Aqua.value.tmp))
    # Type of Sector, i.e. aquaculture or existing, used for the faceting of variables 
    type.name.tmp=x.name;type.name.rep.tmp=x.name
    type.name.rep.tmp[type.name.tmp=="Kelp"|type.name.tmp=="Finfish"|type.name.tmp=="Mussel"]="Aquaculture"
    type.name.rep.tmp[type.name.tmp=="Halibut"|type.name.tmp=="Viewshed"|type.name.tmp=="Benthic"|type.name.tmp=="Disease"]="Impacted Sector"
    Type.of.sector.tmp=rep(type.name.rep.tmp,times=length(Aqua.value.tmp))
    # Type of conventional 
    u.string=rep('Unconstrained',times=length(d.Ux))
    s.string=rep('Unconstrained',times=length(d.Ux))
    Type.of.conventional=c(u.string,s.string)
    value.MSP.out.kelp=data.frame(Aqua.value.tmp,Value.of.MSP.tmp,Group.tmp,Existing.sector.name.tmp,Type.of.sector.tmp,Type.of.conventional)
    Kelp.value.tmp=rbind(Kelp.value.tmp,value.MSP.out.kelp)
  }
}
