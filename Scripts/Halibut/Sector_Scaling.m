if Sector_2==0 && Sector_3==0 && Sector_5==0 && Sector_6==0 &&   Sector_7==0
    %% Do the usual halibut --> write this so it can be used with the 2D (need to consider what is gonna go on if halibut is not used at all) 
else % Other Aquaculture Sectors aka multidimensional --> Work on how to automate this! 
%% Sector 2 
Sector_responseindices=Sector_2indices;
Sector_responseNPVi=Sector_2NPVi;
Sector_2NPV_raw=Sector_2NPVi(Sector_2indices);
X=zeros(size(Sector_indices));
generic_redistrubution_sectors
Sector_2NPVimax=Sector_responseNPVi;
Sector_2NPVmax=Sector_responseNPV;
Sector_responseindices=Sector_2indices;
Sector_responseNPVi=Sector_2NPVi;
X=ones(size(Sector_indices));
generic_redistrubution_sectors
Sector_2NPVimin=Sector_responseNPVi;
Sector_2NPVmin=Sector_responseNPV;
Sector_2NPVi_obj=(Sector_2NPV_raw)./sum(Sector_2NPV_raw);
%% Sector 3 
Sector_responseindices=Sector_3indices;
Sector_responseNPVi=Sector_3NPVi;
Sector_3NPV_raw=Sector_3NPVi(Sector_3indices);
X=zeros(size(Sector_indices));
generic_redistrubution_sectors
Sector_3NPVimax=Sector_responseNPVi;
Sector_3NPVmax=Sector_responseNPV;
Sector_responseindices=Sector_3indices;
Sector_responseNPVi=Sector_3NPVi;
X=ones(size(Sector_indices));
generic_redistrubution_sectors
Sector_3NPVimin=Sector_responseNPVi;
Sector_3NPVmin=Sector_responseNPV;
Sector_3NPVi_obj=(Sector_3NPV_raw)./sum(Sector_3NPV_raw);
%% Sector 4 
Sector_responseindices=Sector_4indices;
Sector_responseNPVi=Sector_4NPVi;
Sector_4NPV_raw=Sector_4NPVi(Sector_4indices);
X=zeros(size(Sector_indices));
generic_redistrubution_sectors
Sector_4NPVimax=Sector_responseNPVi;
Sector_4NPVmax=Sector_responseNPV;
Sector_responseindices=Sector_4indices;
Sector_responseNPVi=Sector_4NPVi;
X=ones(size(Sector_indices));
generic_redistrubution_sectors
Sector_4NPVimin=Sector_responseNPVi;
Sector_4NPVmin=Sector_responseNPV;
Sector_4NPVi_obj=(Sector_4NPV_raw)./sum(Sector_4NPV_raw);
%% Sector 5 
Sector_responseindices=Sector_5indices;
Sector_responseNPVi=Sector_5NPVi;
Sector_5NPV_raw=Sector_5NPVi(Sector_5indices);
X=zeros(size(Sector_indices));
generic_redistrubution_sectors
Sector_5NPVimax=Sector_responseNPVi;
Sector_5NPVmax=Sector_responseNPV;
Sector_responseindices=Sector_5indices;
Sector_responseNPVi=Sector_5NPVi;
X=ones(size(Sector_indices));
generic_redistrubution_sectors
Sector_5NPVimin=Sector_responseNPVi;
Sector_5NPVmin=Sector_responseNPV;
Sector_5NPVi_obj=(Sector_5NPV_raw)./sum(Sector_5NPV_raw);
end