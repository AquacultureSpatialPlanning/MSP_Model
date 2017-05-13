disp('Load data from SeaGrant_allFINALdata_complete_2015.xlsx')
[NUM,TXT,RAW]=xlsread('SeaGrant_allFINALdata_complete_2015.xlsx'); %load spreadsheet with data for all sectors 
col_lat=2; %x-axis = latitude
col_lon=3; %y-axis = longitude
col_FID=1;
target_fid_fulldomain=NUM(:,col_FID);
%% Time dimension of anaylsis 
T=10; %NPV time horizon 
discount_rate=0.05; %discount rate (d=0.05)
discount_factor=1./((1+discount_rate).^[1:T]);
discount_rate_iy_aqua=repmat(discount_factor,length(target_fid_fulldomain),1);
%% Mussel
disp('Parameterizing Mussel...')
disp('Designate Mussel as P1 and V1')
mussel_upfrontcostcol=32; % Load all the mussel variables from the SeaGrant_allFINALdata_complete_2015.xlsx
mussel_yieldkg_col=33;
mussel_annualcost_col=34;
mussel_annualprofit_col=35;
mussel_10yrNPV_est_col=36;
p_mussel=3.3; %Price in dollars per kg
mussel_upfrontcost=NUM(:,mussel_upfrontcostcol); % Cost in the first year per patch
mussel_yieldkg=NUM(:,mussel_yieldkg_col); % Yield per patch per year 
mussel_annualcost=NUM(:,mussel_annualcost_col); % Annual cost per patch for all years following first year 
% mussel_annualprofit=NUM(:,mussel_annualprofit_col); % Annual profit per patch per year (Incorrect rough estimates from Becca, NOT USED IN ANALYSIS) 
% mussel_10yrNPV_est=NUM(:,mussel_10yrNPV_est_col); % Estimated 10 yr NPV per patch  (Incorrect rough estimates from Becca, NOT USED IN ANALYSIS)  
%Calculate true NPV
% Set up the matrices for calculating the total revenue per year per patch 
mussel_rev=mussel_yieldkg.*p_mussel; % Yield of patch multipied by supplied price value 
mussel_rev_mat=repmat(mussel_rev,1,T); % Make a revenue matrix by [# of patches T]
mussel_year_1_total_cost=mussel_upfrontcost+mussel_annualcost; % The first year has to account for the start up costs plus the operating costs for that single year %MAY NEED TO EDIT THIS SINCE THEY MAY NOT PRODUCE ANYTHING IN THE FIRST YEAR 
mussel_cost_mat=[mussel_year_1_total_cost repmat(mussel_annualcost,1,T-1)]; % Make a matrix of all of year 1 cost and a replication of all subsequent years within the range T 
mussel_annual_profit_Piy=mussel_rev_mat-mussel_cost_mat; % Take the difference of the rev matrix and cost matrix to get annual profit per patch per year 
mussel_NPViy=mussel_annual_profit_Piy.*discount_rate_iy_aqua; % Multiply the profit matrix by the discounted rate per year to get the NPV per patch per year 
mussel_NPVi=sum(mussel_NPViy,2); % NPV per patch summed across all years 
mussel_NPV=sum(mussel_NPVi,1); % Total NPV 
% Sorting and scaling 
mussel_NPVi(mussel_NPVi<0)=0; %  Set negative values to zeros 
mussel_NPVi_scaled_soln=mussel_NPVi./sum(mussel_NPVi); % Value scaled for measuring the RELATIVE VALUE OF SOLUTIONS
mussel_NPVi_indices=find(mussel_NPVi>0); % Indices of patches that are postive value 
% P1_raw=mussel_NPVi; % RAW VALUE  
V1=mussel_NPVi>0==1; % Converting zero and postive value to 0 and 1 binary 
P1_full_scaled_obj=mussel_NPVi_scaled_soln; % Scaled for the objective function (Value relative to MIN and MAX impacts)  
%% Fish
disp('Parameterizing FinFish Aquaculture...') 
disp('Designate FinFish as P2 and V2')
Fish_yieldkg_col=41;
Fish_totalcost_col=42;
Fish_annualrev_col=44;
Fish_yieldkg=NUM(:,Fish_yieldkg_col); % Yield per patch per year 
Fish_totalcost=NUM(:,Fish_totalcost_col);  % Annual cost per patch for all years 
Fish_annualrev=NUM(:,Fish_annualrev_col); % Annual revenue per patch for all years 
p_Fish=8; %Price of per kg, NOTE: The set price has a ton to do with what pacthes are set up to be profitable, changing this we estimate will have a dramtic effect on the number and spatial lay out of FISH developable patches, SENSITIVITY ANALYSIS 
% May want to develop a GLM for the price per kg of striped bass 
Fish_rev=Fish_yieldkg.*p_Fish; % Yield of patch multipied by supplied price value 
Fish_rev_mat=repmat(Fish_rev,1,T);  % Make a revenue matrix by [# of patches T]
Fish_cost_mat=[repmat(Fish_totalcost,1,T)]; % Make a replicated matrix with all values per patch of total yearly cost replicated across all years in the value T 
Fish_annual_profit_Piy=Fish_rev_mat-Fish_cost_mat; % Take the difference of the rev matrix and cost matrix to get annual profit per patch per year 
Fish_NPViy=Fish_annual_profit_Piy.*discount_rate_iy_aqua;   % Multiply the profit matrix by the discounted rate per year to get the NPV per patch per year 
Fish_NPVi=sum(Fish_NPViy,2); % NPV per patch summed across all years 
Fish_NPV=sum(Fish_NPVi,1); % Total NPV 
% Sorting and scaling 
Fish_NPVi(Fish_NPVi<0)=0; % Set negative values to zeros 
Fish_NPVi_scaled_soln=Fish_NPVi./sum(Fish_NPVi);  % Value scaled for measuring the RELATIVE VALUE OF SOLUTIONS
Fish_NPVi_indices=find(Fish_NPVi>0); % Indices of patches that are postive value 
% P2_raw=Fish_NPVi; % RAW VALUE  
V2=Fish_NPVi>0==1; % Converting zero and postive value to 0 and 1 binary 
P2_full_scaled_obj=Fish_NPVi_scaled_soln; % Scaled for the objective function (Value relative to MIN and MAX impacts)  
%% Kelp 
disp('loading Kelp...')
disp('Designate Kelp as P3 and V3')
kelp_yieldkg_col=37;
kelp_annualcost_col=38;
kelp_upfrontcost_col=39;
kelp_10yrNPV_est_col=40;
p_kelp=3; %Price in dollars per kg 
kelp_yieldkg=NUM(:,kelp_yieldkg_col);% Yield per patch per year 
kelp_annualcost=NUM(:,kelp_annualcost_col); % Annual cost per patch for all years following first year 
kelp_upfrontcost=NUM(:,kelp_upfrontcost_col);  % Cost in the first year per patch
% kelp_10yrNPV_est=NUM(:,kelp_10yrNPV_est_col);  % Estimated 10 yr NPV per patch  (Incorrect rough estimates from Becca, NOT USED IN ANALYSIS)  
%Calculate true NPV
% Set up the matrices for calculating the total revenue per year per patch 
kelp_rev=kelp_yieldkg.*p_kelp; % Yield of patch multipied by supplied price value 
kelp_rev_mat=repmat(kelp_rev,1,T);  % Make a revenue matrix by [# of patches T]
kelp_year_1_total_cost=kelp_upfrontcost+kelp_annualcost; % The first year has to account for the start up costs plus the operating costs for that single year %MAY NEED TO EDIT THIS SINCE THEY MAY NOT PRODUCE ANYTHING IN THE FIRST YEAR
kelp_cost_mat=[kelp_year_1_total_cost repmat(kelp_annualcost,1,T-1)]; % Make a matrix of all of year 1 cost and a replication of all subsequent years within the range T 
kelp_annual_profit_Piy=kelp_rev_mat-kelp_cost_mat;  % Take the difference of the rev matrix and cost matrix to get annual profit per patch per year 
kelp_NPViy=kelp_annual_profit_Piy.*discount_rate_iy_aqua;   % Multiply the profit matrix by the discounted rate per year to get the NPV per patch per year 
kelp_NPVi=sum(kelp_NPViy,2); % NPV per patch summed across all years 
kelp_NPV=sum(kelp_NPVi,1);  % Total NPV 
% Sorting and scaling 
kelp_NPVi(kelp_NPVi<0)=0;%  Set negative values to zeros 
kelp_NPVi_scaled_soln=kelp_NPVi./sum(kelp_NPVi); % Value scaled for measuring the RELATIVE VALUE OF SOLUTIONS
kelp_NPVi_indices=find(kelp_NPVi>0);  % Indices of patches that are postive value 
% P3_raw=kelp_NPVi;% RAW VALUE  
V3=kelp_NPVi>0==1; % Converting zero and postive value to 0 and 1 binary 
P3_full_scaled_obj=kelp_NPVi_scaled_soln; % Scaled for the objective function (Value relative to MIN and MAX impacts)  
%% Aqua summary parameters 
Aqua_dev_1Y_0N=V1+V2+V3>0==1; % Binary indicating any form of aquaculture 
%% Existing Sectors 
%% Halibut 
disp('loading Halibut...')
disp('Designate Halibut Fishery as P4 and V4')
% Values are most likely outdated (NOT USED IN ANAYLSIS) %NEED TO EDIT FOR THE DUMMY DATA 
Halibut_10yrNPVi=Yi_fulldomain_NPV; %Output from Halibut fishery model 
Halibut_10yrNPVi(Halibut_10yrNPVi<0)=0; %Set negative values to zero 
Halibut_max=sum(Halibut_10yrNPVi);
Halibut_10yrNPV_scaled=Halibut_10yrNPVi./sum(Halibut_10yrNPVi); % Value scaled for measuring the relative value of the solutions
Halibut_10yrNPV_indices=find(Halibut_10yrNPVi>0); % Indices of patches that are postive value
P4_raw=Halibut_10yrNPVi; % RAW VALUE  
V4=Halibut_10yrNPVi>0==1; % Converting zero and postive value to 0 and 1 binary 
P4_full_devi=(Aqua_dev_1Y_0N==0).*P4_raw;% min raw Halibut value, due to full aqua development per patch 
P4_full_dev=sum(P4_full_devi);% min raw Halibut value, due to full aqua development 
P4_aqua_patches=Aqua_dev_1Y_0N.*P4_raw; %All halibut patches that can be impacted by aqua development 
P4_aqua_patches_scaled_obj=P4_aqua_patches./sum(P4_aqua_patches); % Scaled for the objective function (Value relative to MIN and MAX impacts)  
P4_max_dev_scaled_soln=P4_raw./sum(P4_raw); % Scaled for the solutions 

%% Viewshed 
disp('loading Viewshed...');
disp('Designate Viewshed as P5 and V5')
% Create vector of aquadesign 1 is mussel/kelp 2 is finfish 
max_view_impact_aquadesign=zeros(length(target_fid_fulldomain),1);
max_view_impact_aquadesign(V1==1)=1;
max_view_impact_aquadesign(V3==1)=1; % May be overiding mussel cells, does not matter due to assignment of 1 to both 
max_view_impact_aquadesign(V2==1)=2; % Overiding mussel/kelp with finfish where it can go 
Viewshed_impacts_finfish_col_res=46;
Viewshed_impacts_finfish_res=NUM(:,Viewshed_impacts_finfish_col_res); % Residential viewshed impacts for finfish 
Viewshed_impacts_finfish_col_park=49;
Viewshed_impacts_finfish_park=NUM(:,Viewshed_impacts_finfish_col_park); % Recreational viewshed impacts for finfish 
Viewshed_impacts_kelpmussel_col_res=47;
Viewshed_impacts_kelpmussel_res=NUM(:,Viewshed_impacts_kelpmussel_col_res); % Residential viewshed impacts for kelp/mussel
Viewshed_impacts_kelpmussel_col_park=48;
Viewshed_impacts_kelpmussel_park=NUM(:,Viewshed_impacts_kelpmussel_col_park);% Recreational viewshed impacts for kelp/mussel
Viewshed_impacts_finfish=Viewshed_impacts_finfish_res+Viewshed_impacts_finfish_park;% Need to then sum the residential and recreational viewshed columns together for both kind of both kinds of impacted aquaculture 
Viewshed_impacts_kelpmussel=Viewshed_impacts_kelpmussel_res+Viewshed_impacts_kelpmussel_park; %NOTE: Both respective viewshed vectors for mussel/kep and finfish are now the sum of both the residential values and recreational values 
Viewshed_impacts_finfish(Viewshed_impacts_finfish<0)=0;
Viewshed_impacts_kelpmussel(Viewshed_impacts_kelpmussel<0)=0;
Viewshed_max_impact=max([Viewshed_impacts_finfish, Viewshed_impacts_kelpmussel],[],2); % Worst case scenario for viewshed impact 
Viewshed_max_impact_scaled=Viewshed_max_impact/sum(Viewshed_max_impact);
% Scaled the raw impacts wrt the worst case scenario impact 
% Assuming that all patches with postive viewshed impact values can be developed for aquaculture
Viewshed_impacts_finfish_scaled_obj=Viewshed_impacts_finfish./sum(Viewshed_max_impact); % Does not sum to one because mussel/kelp might be possible where finfish isn't 
Viewshed_impacts_kelpmussel_scaled_obj=Viewshed_impacts_kelpmussel./sum(Viewshed_max_impact); % Does not sum to one because even if mussel and kelp development is maxed, impact is not as high as finfish 
Viewshed_impacts_finfish_scaled_soln=Viewshed_impacts_finfish_scaled_obj; % Vector for scaling solutions
Viewshed_impacts_kelpmussel_scaled_soln=Viewshed_impacts_kelpmussel_scaled_obj;% Vector for scaling solutions
% Viewshed_impacts_indices=find(Viewshed_impacts>0);
P5_wrt_finfish=Viewshed_impacts_finfish_scaled_obj; % Scaled for the objective function (Value relative to MIN and MAX impacts)  
P5_wrt_kelpmussel=Viewshed_impacts_kelpmussel_scaled_obj; % Scaled for the objective function (Value relative to MIN and MAX impacts)  
tmp=Viewshed_max_impact>0==1;
V5=tmp;
% tmp1=V1+V2+V3>0==1;
% tmp2=tmp1.*V5;
% P5_full_dev=P5.*tmp2;
% P5_full_dev_scaled_obj=P5_full_dev./sum(P5_full_dev);
%% Enviormental Risk --> This is only affected by Finfish only V6 will be the same everything else will be different based on the values that Becca gives us from Aquamodel 
% disp('loading Enviornmental Risk...');
disp('Designate Enviormental Risk as P6 and V6')
enviro_impact_col=50;
enviro_impact=NUM(:,enviro_impact_col);
fish_environ_logical=zeros(size(Fish_NPVi));
fish_environ_logical(Fish_NPVi_indices)=1; % NOTE: Enviromental Impacts are only associated with patches that are developed for finfish aquaculture 
enviro_impact(enviro_impact<0)=0;
enviro_impact=enviro_impact.*fish_environ_logical;%NOTE: PREVIOUSLY PATCHES THAT HAD ENVIORN IMPACTS WERE INCLUDED EVEN IF THE VALUE OF FISH WAS NEGATIVE!!! 
enviro_impact_scaled_soln=enviro_impact./sum(enviro_impact);
% Assuming that all patches with environmental impacts postive values can be developed for finfish 
P6_aqua_patches_scaled_obj=enviro_impact_scaled_soln;
% P6=enviro_impact_scaled_soln;
%% Aquaculture Indices 
Aqua_dev_indices=find(Aqua_dev_1Y_0N==1);
%% Disease Risk --> Uses the eigenvector centrality matrix to catorgorized each developed patch as hubs or otherwise. "eigenvector centrality takes the entire network into account such that nodes may obtain a high centrality by being connected to many low-centrality nodes or by being connected to a smaller number of high-centrality nodeseigenvector centrality takes the entire network into account such that nodes may obtain a high centrality by being connected to many low-centrality nodes or by being connected to a smaller number of high-centrality nodes"-Randi H. Griffin and Charles L. Nunn  
max_impact_fish_X_logical=V2(Aqua_dev_indices);
max_impact_fish_X=max_impact_fish_X_logical.*3;
X=max_impact_fish_X;
% Calculate Disease Max 
    X_domain=zeros(size(target_fid_fulldomain));
    X_domain(Aqua_dev_indices)=X';
    disease_vector_tmp=[(X_domain==3)+(X_domain==8)+(X_domain==13)+(X_domain==17)];
    disease_vector=disease_vector_tmp(Fish_NPVi_indices);
    tmp1=Fish_NPVi_indices(disease_vector==1); %of the whole domain, cells with aqua dev
    connect_matrix_disease=disease_connect_matrix(tmp1,tmp1);
    adj=connect_matrix_disease;
    eigen_disease=eigencentrality(adj); % Note that eignevalues are sometimes negative, so we need to convert it to postive 
    disease_metrici=abs(eigen_disease);
    disease_metric=sum(disease_metrici);
    disease_maxi=disease_metrici;
    disease_max=disease_metric; 
% Calculate Disease Min 
    eigen_disease=eigencentrality(1); % Feeding a single patch for the minimum, note that eigen centraility will always be 1 for a single node network 
    disease_metrici=abs(eigen_disease);
    disease_metric=sum(disease_metrici);
    disease_mini=disease_metrici;
    disease_min=disease_metric; 
% Plot the maximum disease risk scenario 
X_tmp=X_domain;
X_tmp(X_tmp>0)=1;
X_tmp=logical(X_tmp);
patchmarkersize=14;
patch_color=[1 0 1];
figure
scatter(lat_lon_msp_domain(X_tmp,1),lat_lon_msp_domain(X_tmp,2),patchmarkersize,patch_color,'s','filled') %all patches=grey filled in squares
hold on
scatter(lat_lon_leftmapedge(1),lat_lon_leftmapedge(2),0.1,'.','k')
colorbar
axis tight
xlabel('Latitude')
ylabel('Longitude')
title(['Disease Risk=',num2str(disease_max)])
set(gcf,'color','white'); 
plot(msp_domain(:,1),msp_domain(:,2),'k','linewidth',coastlinewidth) %coastline
box on
axis square
axis tight
%% Parameterizing GA 
%Aqua_dev_indices=find(Aqua_dev_1Y_0N==1);
lb_total=NaN(size(Aqua_dev_indices));
ub_total=NaN(size(Aqua_dev_indices));
% code for each patch 
%[Aqua type possible in patch] [INTCON options in GA] [lb] [ub] 
%      [0,1]                      [0,1]               [0]  [1]
%      [0,2]                      [2,3]               [2]  [3]
%      [0,3]                      [4,5]               [4]  [5]
%      [0,1,2]                    [6,7,8]             [6]  [8]
%      [0,1,3]                    [9,10,11]           [9]  [11]
%      [0,2,3]                    [12,13,14]          [12] [14]
%      [0,1,2,3]                  [15,16,17,18]       [15] [18]

for a=1:length(Aqua_dev_indices)
     if V1(Aqua_dev_indices(a))==1 && V2(Aqua_dev_indices(a))==0 && V3(Aqua_dev_indices(a))==0; lb_total(a)=0; ub_total(a)=1;
    elseif V1(Aqua_dev_indices(a))==0 && V2(Aqua_dev_indices(a))==1 && V3(Aqua_dev_indices(a))==0; lb_total(a)=2; ub_total(a)=3;
    elseif V1(Aqua_dev_indices(a))==0 && V2(Aqua_dev_indices(a))==0 && V3(Aqua_dev_indices(a))==1; lb_total(a)=4; ub_total(a)=5;
    elseif V1(Aqua_dev_indices(a))==1 && V2(Aqua_dev_indices(a))==1 && V3(Aqua_dev_indices(a))==0; lb_total(a)=6; ub_total(a)=8;
    elseif V1(Aqua_dev_indices(a))==1 && V2(Aqua_dev_indices(a))==0 && V3(Aqua_dev_indices(a))==1; lb_total(a)=9; ub_total(a)=11; 
    elseif V1(Aqua_dev_indices(a))==0 && V2(Aqua_dev_indices(a))==1 && V3(Aqua_dev_indices(a))==1; lb_total(a)=12; ub_total(a)=14; 
    elseif V1(Aqua_dev_indices(a))==1 && V2(Aqua_dev_indices(a))==1 && V3(Aqua_dev_indices(a))==1; lb_total(a)=15; ub_total(a)=18; 
    end
end

% Make upper and lower bound for single aquaculture development 
% Mussel 
lb_mussel=NaN(size(Aqua_dev_indices));
ub_mussel=NaN(size(Aqua_dev_indices));
for a=1:length(Aqua_dev_indices)
    if V1(Aqua_dev_indices(a))==1; lb_mussel(a)=0; ub_mussel(a)=1;
    elseif V1(Aqua_dev_indices(a))==0 lb_mussel(a)=0; ub_mussel(a)=0;
    end
end
% Finfish 
lb_finfish=NaN(size(Aqua_dev_indices));
ub_finfish=NaN(size(Aqua_dev_indices));
for a=1:length(Aqua_dev_indices)
    if V2(Aqua_dev_indices(a))==1; lb_finfish(a)=2; ub_finfish(a)=3;
    elseif V2(Aqua_dev_indices(a))==0;  lb_finfish(a)=2; ub_finfish(a)=2;
    end
end
% Kelp 
lb_kelp=NaN(size(Aqua_dev_indices));
ub_kelp=NaN(size(Aqua_dev_indices));
for a=1:length(Aqua_dev_indices)
     if V3(Aqua_dev_indices(a))==1; lb_kelp(a)=4; ub_kelp(a)=5;
     elseif V3(Aqua_dev_indices(a))==0; lb_kelp(a)=4; ub_kelp(a)=4;
    end
end
% Mussel and Finfish 
lb_mussel_finfish=NaN(size(Aqua_dev_indices));
ub_mussel_finfish=NaN(size(Aqua_dev_indices));
for a=1:length(Aqua_dev_indices)
     if V1(Aqua_dev_indices(a))==1 && V2(Aqua_dev_indices(a))==0; lb_mussel_finfish(a)=0; ub_mussel_finfish(a)=1;
    elseif V1(Aqua_dev_indices(a))==0 && V2(Aqua_dev_indices(a))==1; lb_mussel_finfish(a)=2; ub_mussel_finfish(a)=3;
    elseif V1(Aqua_dev_indices(a))==1 && V2(Aqua_dev_indices(a))==1; lb_mussel_finfish(a)=6; ub_mussel_finfish(a)=8;
    elseif V1(Aqua_dev_indices(a))==0 && V2(Aqua_dev_indices(a))==0; lb_mussel_finfish(a)=0; ub_mussel_finfish(a)=0; 
     end
end
% Mussel and Kelp 
lb_mussel_kelp=NaN(size(Aqua_dev_indices));
ub_mussel_kelp=NaN(size(Aqua_dev_indices));
for a=1:length(Aqua_dev_indices)
     if V1(Aqua_dev_indices(a))==1 && V3(Aqua_dev_indices(a))==0; lb_mussel_kelp(a)=0; ub_mussel_kelp(a)=1;
    elseif V1(Aqua_dev_indices(a))==0 && V3(Aqua_dev_indices(a))==1; lb_mussel_kelp(a)=4; ub_mussel_kelp(a)=5;
    elseif V1(Aqua_dev_indices(a))==1 && V3(Aqua_dev_indices(a))==1; lb_mussel_kelp(a)=9; ub_mussel_kelp(a)=11; 
    elseif V1(Aqua_dev_indices(a))==0 && V3(Aqua_dev_indices(a))==0; lb_mussel_kelp(a)=0; ub_mussel_kelp(a)=0; 
    end
end
% Finfish and Kelp 
lb_finfish_kelp=NaN(size(Aqua_dev_indices));
ub_finfish_kelp=NaN(size(Aqua_dev_indices));
for a=1:length(Aqua_dev_indices)
     if V2(Aqua_dev_indices(a))==1 && V3(Aqua_dev_indices(a))==0; lb_finfish_kelp(a)=2; ub_finfish_kelp(a)=3;
    elseif V2(Aqua_dev_indices(a))==0 && V3(Aqua_dev_indices(a))==1; lb_finfish_kelp(a)=4; ub_finfish_kelp(a)=5;
    elseif V2(Aqua_dev_indices(a))==1 && V3(Aqua_dev_indices(a))==1; lb_finfish_kelp(a)=12; ub_finfish_kelp(a)=13; 
    elseif V2(Aqua_dev_indices(a))==0 && V3(Aqua_dev_indices(a))==0; lb_finfish_kelp(a)=0; ub_finfish_kelp(a)=0; 
    end
end

    %% Notes 
    % There is not a single patch in which all three can be in a single [15 16 17 18]
    % patch 
    % Not a single patch which can be just developed for kelp [4 5]
    % Not a single patch which only fish and kelp can be developed [12 13 14]                             
    % Exlude all of these from the if/then statements in the objective
    % function 
%% Quality Control Checks 
% X=zeros(length(Aqua_dev_indices),1);
% [P1,P2,P3,P4,P5,P6,P7,P1i,P2i,P3i,P4i,P5i,P6i,mussel_NPVi_tmp,mussel_NPV_tmp,Fish_NPVi_tmp,Fish_NPV_tmp,kelp_NPVi_tmp,kelp_NPV_tmp,Halibut_NPVi_tmp,Halibut_NPV_tmp,Halibut_NPV_scaled_tmp,Halibut_Yi_tmp,halibut_dynamic_scaled_tmp,Viewshed_raw_valuei_fin_tmp,Viewshed_raw_valuei_kelpmussel_tmp,Viewshed_raw_value_fin_tmp,Viewshed_raw_value_kelpmussel_tmp,enviro_impacti_tmp,enviro_impact_tmp]=Sector_calcs_wrt_X(X,target_fid_fulldomain,Aqua_dev_indices,P1_full_scaled_obj,P2_full_scaled_obj,P3_full_scaled_obj,P4_aqua_patches_scaled_obj,Viewshed_impacts_finfish_scaled_obj,...
%  Viewshed_impacts_kelpmussel_scaled_obj,P6_aqua_patches_scaled_obj,Fish_NPVi_indices,tmp_x,tmp_y,RADII,SQ_1fishable_0notfishable_for_each_soft_depth_patch_ORIGINAL,Nij_initial_msy_TS,x0_PPUE_msy_TS,Fsum_msy,mussel_NPVi,Fish_NPVi,kelp_NPVi,P4_raw,Halibut_max,P4_full_dev,Viewshed_impacts_finfish,Viewshed_impacts_kelpmussel,enviro_impact,target_fid_hab_depth,Wij,numpatches,max_age,T,delta,age_legal,age_mature,Dii,alphaCR,beta_i,habitat_area_i,theta,price,gamma,distance_to_port_for_each_soft_depth_patch,Mii,age_move,Run_full_1_Run_dummy_0,dummy_indices,discount_rate_iy);
% dis
%% Max Impacts 
% General worst case scenario for aquaculture 
    % Calculate effects on each existing sector 
        % Max impacts on halibut 
            Halibut_max_scaled=1;
            Halibut_min_scaled=P4_full_dev/Halibut_max;
            Halibut_value_reduction_at_max_aqua_dev=Halibut_max_scaled-Halibut_min_scaled;
            Halibut_10yrNPVi_tmp=Halibut_10yrNPVi;
            Halibut_10yrNPVi_tmp(Halibut_10yrNPVi_tmp>0)=1;
            Halibut_Y1_N0_no_dev=Halibut_10yrNPVi_tmp; 
            Halibut_Y1_N0_no_dev_logical=logical(Halibut_Y1_N0_no_dev);% Logical matrix showing where halibut is located
            Halibut_Y1_N0_full_dev=Halibut_Y1_N0_no_dev.*(Aqua_dev_1Y_0N==0);
            Halibut_Y1_N0_full_dev_logical=logical(Halibut_Y1_N0_full_dev);% Logical matrix showing where halibut can be fished after full aquaculture development 
        % Max impacts on viewshed 
            max_view_impact_aquadesign_tmp=max_view_impact_aquadesign;
            max_view_impact_aquadesign_tmp(max_view_impact_aquadesign_tmp>0)=1;
            Viewshed_impact_patches=max_view_impact_aquadesign_tmp;
            Viewshed_impact_patches_logical=logical(Viewshed_impact_patches);
            Viewshed_max_impact_sum=sum(Viewshed_max_impact);
            Viewshed_max_impact_scaled_tmp=0;
            Viewshed_min_impact_scaled=1;
        % Max impacts on environment
            enviro_impact_tmp=enviro_impact;
            enviro_impact_tmp(enviro_impact_tmp>0)=1;
            enviro_impact_1Y_0N=enviro_impact_tmp;
            enviro_impact_1Y_0N_logical=logical(enviro_impact_tmp);
            enviro_max_impact_scaled=0;
            enviro_min_impact_scaled=1;
        % Max impacts on disease % Not yet completed 
% Mussel
    % Scenario which yields max impact  
        max_impact_mussel_X=V1(Aqua_dev_indices);
        X=max_impact_mussel_X;
    % Effects on existing sectors 
  [P1,P2,P3,P4,P5,P6,P7,P1i,P2i,P3i,P4i,P5i,P6i,P7i,mussel_NPVi_tmp,mussel_NPV_tmp,Fish_NPVi_tmp,Fish_NPV_tmp,kelp_NPVi_tmp,kelp_NPV_tmp,Halibut_NPVi_tmp,Halibut_NPV_tmp,Halibut_NPV_scaled_tmp,Halibut_Yi_tmp,halibut_dynamic_scaled_tmp,Viewshed_raw_valuei_fin_tmp,Viewshed_raw_valuei_kelpmussel_tmp,Viewshed_raw_value_fin_tmp,Viewshed_raw_value_kelpmussel_tmp,enviro_impacti_tmp,enviro_impact_tmp,disease_rawi,disease_raw_sum]=Sector_calcs_wrt_X(X,target_fid_fulldomain,Aqua_dev_indices,P1_full_scaled_obj,P2_full_scaled_obj,P3_full_scaled_obj,P4_aqua_patches_scaled_obj,Viewshed_impacts_finfish_scaled_obj,...
 Viewshed_impacts_kelpmussel_scaled_obj,P6_aqua_patches_scaled_obj,disease_min,disease_max,Fish_NPVi_indices,SQ_1fishable_0notfishable_for_each_soft_depth_patch_ORIGINAL,Nij_initial_msy_TS,x0_PPUE_msy_TS,Fsum_msy,mussel_NPVi,Fish_NPVi,kelp_NPVi,P4_raw,Halibut_max,P4_full_dev,Viewshed_impacts_finfish,Viewshed_impacts_kelpmussel,enviro_impact,target_fid_hab_depth,Wij,numpatches,max_age,T,delta,age_legal,age_mature,Dii,alphaCR,beta_i,habitat_area_i,theta,price,gamma,distance_to_port_for_each_soft_depth_patch,Mii,age_move,Run_full_1_Run_dummy_0,dummy_indices,discount_rate_iy,disease_connect_matrix);
     % Impacts on Halibut
         Halibut_max_mussel_total_sector_vali=P4i;
         Halibut_max_mussel_total_sector_val=P4;
         Halibut_max_mussel_dynamic_scaled_val=halibut_dynamic_scaled_tmp;
     % Impacts on Viewshed
         Viewshed_max_mussel_total_sector_vali_finfish=Viewshed_raw_valuei_fin_tmp;
         Viewshed_max_mussel_total_sector_vali_mussel_kelp=Viewshed_raw_valuei_kelpmussel_tmp;
         Viewshed_max_mussel_total_sector_val_finfish=Viewshed_raw_value_fin_tmp;
         Viewshed_max_mussel_total_sector_val_mussel_kelp=Viewshed_raw_value_kelpmussel_tmp;
         Viewshed_max_mussel_total_sector_val=P5;
     % Impacts on Enviornmental Impact
         Environmental_max_mussel_total_sector_vali=enviro_impacti_tmp;
         Environmental_max_mussel_total_sector_val=P6;
%          Environmental_max_mussel_total_sector_val=P5;
     % Impacts on Disease 
% Finfish 
    % Scenario which yields max impact  
        max_impact_fish_X_logical=V2(Aqua_dev_indices);
        max_impact_fish_X=max_impact_fish_X_logical.*3;
        X=max_impact_fish_X;
    % Effects on existing sectors 
 [P1,P2,P3,P4,P5,P6,P7,P1i,P2i,P3i,P4i,P5i,P6i,P7i,mussel_NPVi_tmp,mussel_NPV_tmp,Fish_NPVi_tmp,Fish_NPV_tmp,kelp_NPVi_tmp,kelp_NPV_tmp,Halibut_NPVi_tmp,Halibut_NPV_tmp,Halibut_NPV_scaled_tmp,Halibut_Yi_tmp,halibut_dynamic_scaled_tmp,Viewshed_raw_valuei_fin_tmp,Viewshed_raw_valuei_kelpmussel_tmp,Viewshed_raw_value_fin_tmp,Viewshed_raw_value_kelpmussel_tmp,enviro_impacti_tmp,enviro_impact_tmp,disease_rawi,disease_raw_sum]=Sector_calcs_wrt_X(X,target_fid_fulldomain,Aqua_dev_indices,P1_full_scaled_obj,P2_full_scaled_obj,P3_full_scaled_obj,P4_aqua_patches_scaled_obj,Viewshed_impacts_finfish_scaled_obj,...
 Viewshed_impacts_kelpmussel_scaled_obj,P6_aqua_patches_scaled_obj,disease_min,disease_max,Fish_NPVi_indices,SQ_1fishable_0notfishable_for_each_soft_depth_patch_ORIGINAL,Nij_initial_msy_TS,x0_PPUE_msy_TS,Fsum_msy,mussel_NPVi,Fish_NPVi,kelp_NPVi,P4_raw,Halibut_max,P4_full_dev,Viewshed_impacts_finfish,Viewshed_impacts_kelpmussel,enviro_impact,target_fid_hab_depth,Wij,numpatches,max_age,T,delta,age_legal,age_mature,Dii,alphaCR,beta_i,habitat_area_i,theta,price,gamma,distance_to_port_for_each_soft_depth_patch,Mii,age_move,Run_full_1_Run_dummy_0,dummy_indices,discount_rate_iy,disease_connect_matrix);
    % Impacts on Halibut
         Halibut_max_finfish_total_sector_vali=P4i;
         Halibut_max_finfish_total_sector_val=P4;
         Halibut_max_finfish_dynamic_scaled_val=halibut_dynamic_scaled_tmp;
     % Impacts on Viewshed
         Viewshed_max_finfish_total_sector_vali_finfish=Viewshed_raw_valuei_fin_tmp;
         Viewshed_max_finfish_total_sector_vali_mussel_kelp=Viewshed_raw_valuei_kelpmussel_tmp;
         Viewshed_max_finfish_total_sector_val_finfish=Viewshed_raw_value_fin_tmp;
         Viewshed_max_finfish_total_sector_val_mussel_kelp=Viewshed_raw_value_kelpmussel_tmp;
         Viewshed_max_finfish_total_sector_val=P5;
     % Impacts on Enviornmental Impact
         Environmental_max_finfish_total_sector_vali=enviro_impacti_tmp;
         Environmental_max_finfish_total_sector_val=P6;
     % Impacts on Disease 
        Disease_max_finfish_total_sector_vali=P7i;
         Disease_max_finfish_total_sector_val=P7;
% Kelp 
    % Scenario which yields max impact  
        max_impact_kelp_X_logical=V3(Aqua_dev_indices);
        max_impact_kelp_X=max_impact_kelp_X_logical.*5;
        X=max_impact_kelp_X;
    % Effects on existing sectors 
[P1,P2,P3,P4,P5,P6,P7,P1i,P2i,P3i,P4i,P5i,P6i,P7i,mussel_NPVi_tmp,mussel_NPV_tmp,Fish_NPVi_tmp,Fish_NPV_tmp,kelp_NPVi_tmp,kelp_NPV_tmp,Halibut_NPVi_tmp,Halibut_NPV_tmp,Halibut_NPV_scaled_tmp,Halibut_Yi_tmp,halibut_dynamic_scaled_tmp,Viewshed_raw_valuei_fin_tmp,Viewshed_raw_valuei_kelpmussel_tmp,Viewshed_raw_value_fin_tmp,Viewshed_raw_value_kelpmussel_tmp,enviro_impacti_tmp,enviro_impact_tmp,disease_rawi,disease_raw_sum]=Sector_calcs_wrt_X(X,target_fid_fulldomain,Aqua_dev_indices,P1_full_scaled_obj,P2_full_scaled_obj,P3_full_scaled_obj,P4_aqua_patches_scaled_obj,Viewshed_impacts_finfish_scaled_obj,...
 Viewshed_impacts_kelpmussel_scaled_obj,P6_aqua_patches_scaled_obj,disease_min,disease_max,Fish_NPVi_indices,SQ_1fishable_0notfishable_for_each_soft_depth_patch_ORIGINAL,Nij_initial_msy_TS,x0_PPUE_msy_TS,Fsum_msy,mussel_NPVi,Fish_NPVi,kelp_NPVi,P4_raw,Halibut_max,P4_full_dev,Viewshed_impacts_finfish,Viewshed_impacts_kelpmussel,enviro_impact,target_fid_hab_depth,Wij,numpatches,max_age,T,delta,age_legal,age_mature,Dii,alphaCR,beta_i,habitat_area_i,theta,price,gamma,distance_to_port_for_each_soft_depth_patch,Mii,age_move,Run_full_1_Run_dummy_0,dummy_indices,discount_rate_iy,disease_connect_matrix);
     % Impacts on Halibut 
     % Impacts on Halibut
         Halibut_max_kelp_total_sector_vali=P4i;
         Halibut_max_kelp_total_sector_val=P4;
         Halibut_max_kelp_dynamic_scaled_val=halibut_dynamic_scaled_tmp;
     % Impacts on Viewshed
         Viewshed_max_kelp_total_sector_vali_finfish=Viewshed_raw_valuei_fin_tmp;
         Viewshed_max_kelp_total_sector_vali_mussel_kelp=Viewshed_raw_valuei_kelpmussel_tmp;
         Viewshed_max_kelp_total_sector_val_finfish=Viewshed_raw_value_fin_tmp;
         Viewshed_max_kelp_total_sector_val_mussel_kelp=Viewshed_raw_value_kelpmussel_tmp;
         Viewshed_max_kelp_total_sector_val=P5;
     % Impacts on Enviornmental Impact
         Environmental_max_kelp_total_sector_vali=enviro_impacti_tmp;
         Environmental_max_kelp_total_sector_val=P6;
% Mussel and Finfish
% Scenario which yields max impact  
        max_impact_mussel_X=V1(Aqua_dev_indices);
        max_impact_mussel_and_finfish_tmp1=max_impact_mussel_X.*7;
        max_impact_fish_X_logical=V2(Aqua_dev_indices);
        max_impact_mussel_and_finfish_tmp2=max_impact_fish_X_logical.*8;
        max_impact_mussel_and_finfish_tmp=max_impact_mussel_and_finfish_tmp1+max_impact_mussel_and_finfish_tmp2;
        max_impact_mussel_and_finfish_tmp3(max_impact_mussel_and_finfish_tmp==15)=7;
        max_impact_mussel_and_finfish_heavy_mussel=max_impact_mussel_and_finfish_tmp3;
        max_impact_mussel_and_finfish_tmp4(max_impact_mussel_and_finfish_tmp==15)=8;
        max_impact_mussel_and_finfish_heavy_finfish=max_impact_mussel_and_finfish_tmp4;
% Mussel and Kelp 
        max_impact_mussel_X=V1(Aqua_dev_indices);
        max_impact_mussel_and_kelp_tmp1=max_impact_mussel_X.*10;
        max_impact_kelp_X_logical=V3(Aqua_dev_indices);
        max_impact_mussel_and_kelp_tmp2=max_impact_kelp_X_logical.*11;
        max_impact_mussel_and_kelp_tmp=max_impact_mussel_and_kelp_tmp1+max_impact_mussel_and_kelp_tmp2;
        max_impact_mussel_and_kelp_tmp3(max_impact_mussel_and_kelp_tmp==21)=10;
        max_impact_mussel_and_kelp_heavy_mussel=max_impact_mussel_and_kelp_tmp3;
        max_impact_mussel_and_kelp_tmp4(max_impact_mussel_and_kelp_tmp==21)=11;
        max_impact_mussel_and_kelp_heavy_kelp=max_impact_mussel_and_kelp_tmp4;
%         X=max_impact_mussel_X;
%         max_impact_mussel_and_finfish_X_logical=V2(Aqua_dev_indices);
%         max_impact_fish_X=max_impact_fish_X_logical.*3;
%         X=max_impact_fish_X;
    % Effects on existing sectors 
%  [P1,P2,P3,P4,P5,P6,P7,P1i,P2i,P3i,P4i,P5i,P6i,P7i,mussel_NPVi_tmp,mussel_NPV_tmp,Fish_NPVi_tmp,Fish_NPV_tmp,kelp_NPVi_tmp,kelp_NPV_tmp,Halibut_NPVi_tmp,Halibut_NPV_tmp,Halibut_NPV_scaled_tmp,Halibut_Yi_tmp,halibut_dynamic_scaled_tmp,Viewshed_raw_valuei_fin_tmp,Viewshed_raw_valuei_kelpmussel_tmp,Viewshed_raw_value_fin_tmp,Viewshed_raw_value_kelpmussel_tmp,enviro_impacti_tmp,enviro_impact_tmp,disease_rawi,disease_raw_sum]=Sector_calcs_wrt_X(X,target_fid_fulldomain,Aqua_dev_indices,P1_full_scaled_obj,P2_full_scaled_obj,P3_full_scaled_obj,P4_aqua_patches_scaled_obj,Viewshed_impacts_finfish_scaled_obj,...
%  Viewshed_impacts_kelpmussel_scaled_obj,P6_aqua_patches_scaled_obj,disease_min,disease_max,Fish_NPVi_indices,SQ_1fishable_0notfishable_for_each_soft_depth_patch_ORIGINAL,Nij_initial_msy_TS,x0_PPUE_msy_TS,Fsum_msy,mussel_NPVi,Fish_NPVi,kelp_NPVi,P4_raw,Halibut_max,P4_full_dev,Viewshed_impacts_finfish,Viewshed_impacts_kelpmussel,enviro_impact,target_fid_hab_depth,Wij,numpatches,max_age,T,delta,age_legal,age_mature,Dii,alphaCR,beta_i,habitat_area_i,theta,price,gamma,distance_to_port_for_each_soft_depth_patch,Mii,age_move,Run_full_1_Run_dummy_0,dummy_indices,discount_rate_iy,disease_connect_matrix);
%     
%% Plot results 
    % Map of patches that can be developed for aquaculture 
    patchmarkersize=14;
            patch_color=[1 0 1];
            figure
            scatter(lat_lon_msp_domain(Aqua_dev_1Y_0N,1),lat_lon_msp_domain(Aqua_dev_1Y_0N,2),patchmarkersize,patch_color,'s','filled') %all patches=grey filled in squares
            hold on
            scatter(lat_lon_leftmapedge(1),lat_lon_leftmapedge(2),0.1,'.','k')
            % colorbar
            axis tight
            xlabel('Latitude')
            ylabel('Longitude')
            title('Developable Patches for Aquaculture')
            set(gcf,'color','white'); 
            plot(msp_domain(:,1),msp_domain(:,2),'k','linewidth',coastlinewidth) %coastline
            box on
            % axis square
            axis tight
            savefig('Aqua_patches')
    % Map of patches that can be developed for specific kind of aqua 
        % Mussel 
            patch_color=[0 0 1];
            figure
            scatter(lat_lon_msp_domain(V1,1),lat_lon_msp_domain(V1,2),patchmarkersize,patch_color,'s','filled') %all patches=grey filled in squares
            hold on
            scatter(lat_lon_leftmapedge(1),lat_lon_leftmapedge(2),0.1,'.','k')
            % colorbar
            axis tight
            xlabel('Latitude')
            ylabel('Longitude')
            title('Developable Patches for Mussel')
            set(gcf,'color','white'); 
            plot(msp_domain(:,1),msp_domain(:,2),'k','linewidth',coastlinewidth) %coastline
            box on
            % axis square
            axis tight
            savefig('Mussel_patches')
            
            patch_color=mussel_NPVi_scaled_soln(V1);
            figure
            tmp=lat_lon_msp_domain(V1,:);
            scatter(tmp(:,1),tmp(:,2),patchmarkersize,patch_color,'s','filled') %all patches=grey filled in squares
            hold on
            scatter(lat_lon_leftmapedge(1),lat_lon_leftmapedge(2),0.1,'.','k')
            % colorbar
            axis tight
            xlabel('Latitude')
            ylabel('Longitude')
            title('Developable Patches for Mussel')
            set(gcf,'color','white'); 
            plot(msp_domain(:,1),msp_domain(:,2),'k','linewidth',coastlinewidth) %coastline
            box on
            % axis square
            axis tight
            savefig('Mussel_patches_w_value')
        % Finfish 
            patch_color=[1 0 0];
            figure
            scatter(lat_lon_msp_domain(V2,1),lat_lon_msp_domain(V2,2),patchmarkersize,patch_color,'s','filled') %all patches=grey filled in squares
            hold on
            scatter(lat_lon_leftmapedge(1),lat_lon_leftmapedge(2),0.1,'.','k')
            % colorbar
            axis tight
            xlabel('Latitude')
            ylabel('Longitude')
            title('Developable Patches for Finfish')
            set(gcf,'color','white'); 
            plot(msp_domain(:,1),msp_domain(:,2),'k','linewidth',coastlinewidth) %coastline
            box on
            % axis square
            axis tight
            savefig('Finfish_patches')
            
            patch_color=Fish_NPVi_scaled_soln(V2);
            tmp=lat_lon_msp_domain(V2,:);
            figure
            scatter(tmp(:,1),tmp(:,2),patchmarkersize,patch_color,'s','filled') %all patches=grey filled in squares
            hold on
            scatter(lat_lon_leftmapedge(1),lat_lon_leftmapedge(2),0.1,'.','k')
            % colorbar
            axis tight
            xlabel('Latitude')
            ylabel('Longitude')
            title('Developable Patches for Finfish')
            set(gcf,'color','white'); 
            plot(msp_domain(:,1),msp_domain(:,2),'k','linewidth',coastlinewidth) %coastline
            box on
            % axis square
            axis tight
            savefig('Finfish_patches_w_value')
        % Kelp 
            patch_color=[0 1 0];
            figure
            scatter(lat_lon_msp_domain(V3,1),lat_lon_msp_domain(V3,2),patchmarkersize,patch_color,'s','filled') %all patches=grey filled in squares
            hold on
            scatter(lat_lon_leftmapedge(1),lat_lon_leftmapedge(2),0.1,'.','k')
            % colorbar
            axis tight
            xlabel('Latitude')
            ylabel('Longitude')
            title('Developable Patches for Kelp')
            set(gcf,'color','white'); 
            plot(msp_domain(:,1),msp_domain(:,2),'k','linewidth',coastlinewidth) %coastline
            box on
            % axis square
            axis tight
            savefig('Kelp_patches')
            
            patch_color=kelp_NPVi_scaled_soln(V3);
            tmp=lat_lon_msp_domain(V3,:);
            figure
            scatter(tmp(:,1),tmp(:,2),patchmarkersize,patch_color,'s','filled') %all patches=grey filled in squares
            hold on
            scatter(lat_lon_leftmapedge(1),lat_lon_leftmapedge(2),0.1,'.','k')
            % colorbar
            axis tight
            xlabel('Latitude')
            ylabel('Longitude')
            title('Developable Patches for Kelp')
            set(gcf,'color','white'); 
            plot(msp_domain(:,1),msp_domain(:,2),'k','linewidth',coastlinewidth) %coastline
            box on
            % axis square
            axis tight
            savefig('Kelp_patches_w_value')
    % Map of existing sector area that can be impacted by aquaculture 
        % Halibut 
            %Plot of aquaculture developable patches
            patch_color=[1 0 1];
            figure
            scatter(lat_lon_msp_domain(Halibut_Y1_N0_no_dev_logical,1),lat_lon_msp_domain(Halibut_Y1_N0_no_dev_logical,2),patchmarkersize,patch_color,'s','filled') %all patches=grey filled in squares
            hold on
            scatter(lat_lon_leftmapedge(1),lat_lon_leftmapedge(2),0.1,'.','k')
            % colorbar
            axis tight
            xlabel('Latitude')
            ylabel('Longitude')
            title('Halibut Patches No Aqua Development')
            set(gcf,'color','white'); 
            plot(msp_domain(:,1),msp_domain(:,2),'k','linewidth',coastlinewidth) %coastline
            box on
            % axis square
            axis tight
            savefig('Halibut_no_dev_map')
            
            patch_color=[1 0 1];
            figure
            scatter(lat_lon_msp_domain(Halibut_Y1_N0_full_dev_logical,1),lat_lon_msp_domain(Halibut_Y1_N0_full_dev_logical,2),patchmarkersize,patch_color,'s','filled') %all patches=grey filled in squares
            hold on
            scatter(lat_lon_leftmapedge(1),lat_lon_leftmapedge(2),0.1,'.','k')
            % colorbar
            axis tight
            xlabel('Latitude')
            ylabel('Longitude')
            title('Halibut Patches Full Aqua Development')
            set(gcf,'color','white'); 
            plot(msp_domain(:,1),msp_domain(:,2),'k','linewidth',coastlinewidth) %coastline
            box on
            % axis square
            axis tight
            savefig('Halibut_full_dev_map')
            
            patch_color=P4_max_dev_scaled_soln(V4);
            tmp=lat_lon_msp_domain(V4,:);
            figure
            scatter(tmp(:,1),tmp(:,2),patchmarkersize,patch_color,'s','filled') %all patches=grey filled in squares
            hold on
            scatter(lat_lon_leftmapedge(1),lat_lon_leftmapedge(2),0.1,'.','k')
            % colorbar
            axis tight
            xlabel('Latitude')
            ylabel('Longitude')
            title('Economic Value of Halibut Patches')
            set(gcf,'color','white'); 
            plot(msp_domain(:,1),msp_domain(:,2),'k','linewidth',coastlinewidth) %coastline
            box on
            % axis square
            axis tight
            savefig('Halibut_value_map')
            
            patch_color=P4_max_dev_scaled_soln(Halibut_Y1_N0_full_dev_logical);
            tmp=lat_lon_msp_domain(Halibut_Y1_N0_full_dev_logical,:);
            figure
            scatter(tmp(:,1),tmp(:,2),patchmarkersize,patch_color,'s','filled') %all patches=grey filled in squares
            hold on
            scatter(lat_lon_leftmapedge(1),lat_lon_leftmapedge(2),0.1,'.','k')
            % colorbar
            axis tight
            xlabel('Latitude')
            ylabel('Longitude')
            title('Economic Value of Halibut Patches')
            set(gcf,'color','white'); 
            plot(msp_domain(:,1),msp_domain(:,2),'k','linewidth',coastlinewidth) %coastline
            box on
            % axis square
            axis tight
            savefig('Halibut_value_map_with_development')
        % Viewshed 
            patch_color=[0 1 1];
            figure
            scatter(lat_lon_msp_domain(Viewshed_impact_patches_logical,1),lat_lon_msp_domain(Viewshed_impact_patches_logical,2),patchmarkersize,patch_color,'s','filled') %all patches=grey filled in squares
            hold on
            scatter(lat_lon_leftmapedge(1),lat_lon_leftmapedge(2),0.1,'.','k')
            % colorbar
            axis tight
            xlabel('Latitude')
            ylabel('Longitude')
            title('Viewshed Impact Patches at Full Aqua Development')
            set(gcf,'color','white'); 
            plot(msp_domain(:,1),msp_domain(:,2),'k','linewidth',coastlinewidth) %coastline
            box on
            % axis square
            axis tight
            savefig('Viewshed_impacts_map')
            
            patch_color=Viewshed_max_impact_scaled(Aqua_dev_1Y_0N);
            tmp=lat_lon_msp_domain(Aqua_dev_1Y_0N,:);
            figure
            scatter(tmp(:,1),tmp(:,2),patchmarkersize,patch_color,'s','filled') %all patches=grey filled in squares
            hold on
            scatter(lat_lon_leftmapedge(1),lat_lon_leftmapedge(2),0.1,'.','k')
            % colorbar
            axis tight
            xlabel('Latitude')
            ylabel('Longitude')
            title('Viewshed Impact Patches at Full Aqua Development')
            set(gcf,'color','white'); 
            plot(msp_domain(:,1),msp_domain(:,2),'k','linewidth',coastlinewidth) %coastline
            box on
            % axis square
            axis tight
            savefig('Viewshed_impacts_map')
        % Enviornmental Impacts 
            patch_color=[1 0 0];
            figure
            scatter(lat_lon_msp_domain(enviro_impact_1Y_0N_logical,1),lat_lon_msp_domain(enviro_impact_1Y_0N_logical,2),patchmarkersize,patch_color,'s','filled') %all patches=grey filled in squares
            hold on
            scatter(lat_lon_leftmapedge(1),lat_lon_leftmapedge(2),0.1,'.','k')
            % colorbar
            axis tight
            xlabel('Latitude')
            ylabel('Longitude')
            title('Environmental Impact Patches at Full Aqua Development')
            set(gcf,'color','white'); 
            plot(msp_domain(:,1),msp_domain(:,2),'k','linewidth',coastlinewidth) %coastline
            box on
            % axis square
            axis tight
            savefig('Environmental_impacts_map')
            
            patch_color=enviro_impact_scaled_soln(V2);
            tmp=lat_lon_msp_domain(V2,:);
            figure
            scatter(tmp(:,1),tmp(:,2),patchmarkersize,patch_color,'s','filled') %all patches=grey filled in squares
            hold on
            scatter(lat_lon_leftmapedge(1),lat_lon_leftmapedge(2),0.1,'.','k')
            % colorbar
            axis tight
            xlabel('Latitude')
            ylabel('Longitude')
            title('Environmental Impact Patches at Full Aqua Development')
            set(gcf,'color','white'); 
            plot(msp_domain(:,1),msp_domain(:,2),'k','linewidth',coastlinewidth) %coastline
            box on
            % axis square
            axis tight
            savefig('Environmental_impacts_map_w_value')
    % Graph
close all 





