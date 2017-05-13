function [P1,P2,P3,P4,P5,P6,P7,P1i,P2i,P3i,P4i,P5i,P6i,P7i,mussel_NPVi_tmp,mussel_NPV_tmp,Fish_NPVi_tmp,Fish_NPV_tmp,kelp_NPVi_tmp,kelp_NPV_tmp,Halibut_NPVi_tmp,Halibut_NPV_tmp,Halibut_NPV_scaled_tmp,Halibut_Yi_tmp,halibut_dynamic_scaled_tmp,Viewshed_raw_valuei_fin_tmp,Viewshed_raw_valuei_kelpmussel_tmp,Viewshed_raw_value_fin_tmp,Viewshed_raw_value_kelpmussel_tmp,enviro_impacti_tmp,enviro_impact_tmp,disease_rawi,disease_raw_sum]=Sector_calcs_wrt_X(X,target_fid_fulldomain,Aqua_dev_indices,P1_full_scaled_obj,P2_full_scaled_obj,P3_full_scaled_obj,P4_aqua_patches_scaled_obj,Viewshed_impacts_finfish_scaled_obj,...
 Viewshed_impacts_kelpmussel_scaled_obj,P6_aqua_patches_scaled_obj,disease_min,disease_max,Fish_NPVi_indices,SQ_1fishable_0notfishable_for_each_soft_depth_patch_ORIGINAL,Nij_initial_msy_TS,x0_PPUE_msy_TS,Fsum_msy,mussel_NPVi,Fish_NPVi,kelp_NPVi,P4_raw,Halibut_max,P4_full_dev,Viewshed_impacts_finfish,Viewshed_impacts_kelpmussel,enviro_impact,target_fid_hab_depth,Wij,numpatches,max_age,T,delta,age_legal,age_mature,Dii,alphaCR,beta_i,habitat_area_i,theta,price,gamma,distance_to_port_for_each_soft_depth_patch,Mii,age_move,Run_full_1_Run_dummy_0,dummy_indices,discount_rate_iy,disease_connect_matrix)
X_domain=zeros(size(target_fid_fulldomain));
 X_domain(Aqua_dev_indices)=X';
 % These are the patch specific sector values 

P1i=P1_full_scaled_obj.*(X_domain==1);
    mussel_NPVi_tmp=mussel_NPVi.*(P1i>0==1);
    mussel_NPV_tmp=sum(mussel_NPVi(P1i>0));
P2i=P2_full_scaled_obj.*(X_domain==2);
    Fish_NPVi_tmp=Fish_NPVi.*(P2i>0==1);
    Fish_NPV_tmp=sum(Fish_NPVi(P2i>0));
P3i=P3_full_scaled_obj.*(X_domain==3);
    kelp_NPVi_tmp=kelp_NPVi.*(P3i>0==1);
    kelp_NPV_tmp=sum(kelp_NPVi(P3i>0));
P4i=P4_aqua_patches_scaled_obj.*(X_domain==0);
    Halibut_NPVi_tmp=P4_raw.*(P1i+P2i+P3i>0==0);
    Halibut_NPV_tmp=sum(P4_raw(P1i+P2i+P3i==0));
    Halibut_NPV_scaled_tmp=(Halibut_NPV_tmp-P4_full_dev)/(Halibut_max-P4_full_dev); % May want to edit this! 
P5i=(Viewshed_impacts_finfish_scaled_obj.*(X_domain==2))+(Viewshed_impacts_kelpmussel_scaled_obj.*(X_domain==1)+(X_domain==3));% Need to edit 
    Viewshed_raw_valuei_fin_tmp=Viewshed_impacts_finfish.*(P2i>0==1);
    Viewshed_raw_valuei_kelpmussel_tmp=Viewshed_impacts_kelpmussel.*(P1i+P3i>0==1);
    Viewshed_raw_value_fin_tmp=sum(Viewshed_impacts_finfish(P2i>0));
    Viewshed_raw_value_kelpmussel_tmp=sum(Viewshed_impacts_kelpmussel(P1i+P3i>0));
P6i=P6_aqua_patches_scaled_obj.*(X_domain==2);
    enviro_impacti_tmp=enviro_impact.*(P6i>0==1);
    enviro_impact_tmp=sum(enviro_impact(P6i>0));

% Total Sector values 
P1=sum(P1i);
P2=sum(P2i);
P3=sum(P3i);
P4=sum(P4i);
P5_low=sum(P5i);
P6_low=sum(P6i);

%% Halibut 
SQ_1fishable_0notfishable_for_each_soft_depth_patch=SQ_1fishable_0notfishable_for_each_soft_depth_patch_ORIGINAL; %Reset fishable vector
tmp1=P1i+P2i+P3i>0; %find which aqua patches were developed
tmp2=target_fid_fulldomain(tmp1); %figure out what fid numbers in the full domain represent those patches
[C, IA, IB]=intersect(tmp2,target_fid_hab_depth); %find fid's in common with those that represent the halibut patches. IB=indices (not fid's) of the halibut patches that are developed by aqua dev
SQ_1fishable_0notfishable_for_each_soft_depth_patch(IB)=0; %set those patches to not fishable (note: they may already be not fishable due to other constraints)
%Run TS model
Nij_initial=Nij_initial_msy_TS; %initial stock status (patch and age class specific)
x0_PPUE_initial=x0_PPUE_msy_TS;
Fsum=Fsum_msy; %TAE (fleet model) <-Assuming no change in management, only in fishery distribution. 
%I think the above assumption is safe because at most aqua closes
%down X% of halibut fishing grounds (what is X?), thus optimal TAE will probably not be much different
%Determine dynamic fishery result using fleet model or Ei displacement model
Ei_fleet1_displace2=1;
if Ei_fleet1_displace2==1
    BasicSpatialHarvestModel_COREv1_fleet_TS
elseif Ei_fleet1_displace2==1
    Ei=Ei_MPAdisplace_proportional_to_fishable_patches(SQ_1fishable_0notfishable_for_each_soft_depth_patch,Ei_wrt_Fsum_msy);
    Eiy=repmat(Ei,1,y);
    BasicSpatialHarvestModel_COREv1_TS 
end
%Calc and record results
%calc patch specific NPV of fishery
% if Run_full_1_Run_dummy_0==0
%     tmp1=Yiy(dummy_indices,:).*discount_rate_iy(dummy_indices,:); %multiply year and patch specific values by discounted value
% else
tmp1=Yiy.*discount_rate_iy; %discount_rate_iy_halibut_Ei
% end
Yi_NPV_wrt_alpha=sum(tmp1,2); %sum years together to calc NPV for each patch
Halibut_Yi_tmp=sum(Yi_NPV_wrt_alpha);
halibut_dynamic_scaled_tmp=(Halibut_Yi_tmp-P4_full_dev)/(Halibut_max-P4_full_dev);


%% Disease 
% disease_vector_tmp=[(X_domain==3)+(X_domain==8)+(X_domain==13)+(X_domain==17)]; % Integers which designate fish aqua development i.e. a disease risk 
% disease_vector=disease_vector_tmp(Fish_NPVi_indices);
% % Network Disease Approach 
% tmp1=Fish_NPVi_indices(disease_vector==1); %of the fish aquaculture domain, cells with are developed given the spatial plan 
% connect_matrix_disease=disease_connect_matrix(tmp1,tmp1);
% adj=connect_matrix_disease;
% eigen_disease=eigencentrality(adj); % Note that eignevalues are always negative, so we need to convert it to postive 
% disease_rawi=abs(eigen_disease); 
% disease_raw_sum=sum(disease_rawi);
% %P7i=disease_rawi./disease_raw_sum; %How each individual node influences the total eigenvector centrality sum 
% disease_raw_sum(disease_raw_sum==0)=disease_min;
% P7_low=(disease_raw_sum-disease_min)/(disease_max-disease_min);
% P7i_dev=(disease_rawi./disease_raw_sum).*P7_low; %How each individual node influences the total eigenvector centrality sum  
% P7i_dev(isempty(P7i_dev))=0;
% tmp1=Fish_NPVi_indices(P7i_dev>0); %of the whole domain, cells with aqua dev
% tmp2=Fish_NPVi_indices;
% [C, IA, IB]=intersect(tmp1,tmp2);
% P7i=zeros(size(target_fid_fulldomain,1),1);
% P7i(C)=P7i_dev;
disease_vector_tmp=(X_domain==2); % Integers which designate fish aqua development i.e. a disease risk 
disease_vector=disease_vector_tmp(Fish_NPVi_indices);
% Network Disease Approach 
tmp1=Fish_NPVi_indices(disease_vector==1); %of the fish aquaculture domain, cells with are developed given the spatial plan 
connect_matrix_disease=disease_connect_matrix(tmp1,tmp1);
adj=connect_matrix_disease;
eigen_disease=eigencentrality(adj); % Note that eignevalues are always negative, so we need to convert it to postive 
disease_rawi=abs(eigen_disease); 
disease_raw_sum=sum(disease_rawi);
% P7i=disease_rawi./disease_raw_sum; %How each individual node influences the total eigenvector centrality sum 
disease_raw_sum(disease_raw_sum==0)=disease_min;
P7_low=(disease_raw_sum-disease_min)/(disease_max-disease_min);
P7i_dev=(disease_rawi./disease_raw_sum).*P7_low; %How each individual node influences the total eigenvector centrality sum  
P7i_dev(isempty(P7i_dev))=0;
if sum(P7i_dev)==0
    P7i=zeros(size(target_fid_fulldomain,1),1);
else
% tmp1=Fish_NPVi_indices(P7i_dev>0); %of the whole domain, cells with aqua dev
% tmp2=Fish_NPVi_indices;
% [C, IA, IB]=intersect(tmp1,tmp2);
P7i=zeros(size(target_fid_fulldomain,1),1);
P7i(tmp1)=P7i_dev;
end% P7i(C)=P7i_dev;
%P7i=P7i_dev; %NOTE: NEED TO MAKE THIS THE SAME LENGTH AS OTHER SECTORS Pni
%P7i=disease_metrici;
%else 
%disease_metric=0;
%P7i=0;
%end
% disease_riski_tmp=disease_metrici;
% disease_risk_tmp=disease_metric;
% P7i=P7_aqua_patches_scaled_obj.*disease_vector;
% P7_low=sum(P7i);
% Transform the sectors where high values are good: Viewshed, Enviormental Impacts and Disease  
P5=abs(1-P5_low);
P6=abs(1-P6_low);
P7=abs(1-P7_low);