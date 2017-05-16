%Calcuates the spatial pattern of fishing effort (Ei) given PPUE=Average profit = profit per unit effort
function [out]=Ei_wrt_PPUE(x0_PPUE,habitat_area_i,theta,price,Bi_legal,delta,gamma,distance_to_port_for_each_soft_depth_patch,SQ_1fishable_0notfishable_for_each_soft_depth_patch)

% whos x0_PPUE habitat_area_i theta price Bi_legal delta gamma distance_to_port_for_each_soft_depth_patch SQ_1fishable_0notfishable_for_each_soft_depth_patch

% global Bi_legal price delta theta Fsum habitat_area_i gamma distance_to_port_for_each_soft_depth_patch SQ_1fishable_0notfishable_for_each_soft_depth_patch qi
 load('Tuner_save_Jun_22_2015_15_28','habitat_area_i','theta','price','Bi_legal','delta', 'gamma', 'distance_to_port_for_each_soft_depth_patch')                   

%      %(EXP MODEL) w/o distance from port effect
%      Z=(x0_PPUE+habitat_area_i.*theta)./(price.*Bi_legal);
%     lambert_z=-exp(-1./Z)./Z;
%     %wapr will break if lambert_z is out of range (i.e., <-exp(-1)); so set
%     %those values to the edge of the range = -exp(-1)
%     lambert_z(lambert_z<-exp(-1))=-exp(-1);
%     [lambert_w nerror]=wapr(lambert_z); %use stock Matlab function lambertw once get newer version
%     deltai_plus_Ei=1./Z+lambert_w';
%     deltai_plus_Ei=max([zeros(length(deltai_plus_Ei),1) deltai_plus_Ei],[],2); %get ride of tiny negative results
%     Ei_vector_fleet=deltai_plus_Ei-delta; %remove natural mortality effect
%     Ei_vector_fleet=max(Ei_vector_fleet,0); %only count positive harvest patches
%     out=Ei_vector_fleet;

%      %(EXP MODEL) WITH distance from port effect
%      Z=(x0_PPUE+habitat_area_i.*theta.*(1+gamma.*distance_to_port_for_each_soft_depth_patch))./(price.*Bi_legal);
%     lambert_z=-exp(-1./Z)./Z;
%     %wapr will break if lambert_z is out of range (i.e., <-exp(-1)); so set
%     %those values to the edge of the range = -exp(-1)
%     lambert_z(lambert_z<-exp(-1))=-exp(-1);
%     [lambert_w nerror]=wapr(lambert_z); %use stock Matlab function lambertw once get newer version
%     deltai_plus_Ei=1./Z+lambert_w';
%     deltai_plus_Ei=max([zeros(length(deltai_plus_Ei),1) deltai_plus_Ei],[],2); %get ride of tiny negative results
%     Ei_vector_fleet=deltai_plus_Ei-delta; %remove natural mortality effect
%     Ei_vector_fleet=max(Ei_vector_fleet,0); %only count positive harvest patches
%     out=Ei_vector_fleet;

    %ACCOUNTING FOR SQ FISH EXCLUSION PATHCES
     %(EXP MODEL) WITH distance from port effect
     Ei_vector_fleet_all_patches=zeros(size(habitat_area_i)); %This will be the output of Ei for ALL patches, fill in with fleet model results from FISHABLE patches
     %pull just the fishable patches:
     Bi_legal_fishable=Bi_legal(SQ_1fishable_0notfishable_for_each_soft_depth_patch); 
     habitat_area_i_fishable=habitat_area_i(SQ_1fishable_0notfishable_for_each_soft_depth_patch);
     distance_to_port_for_each_soft_depth_patch_fishable=distance_to_port_for_each_soft_depth_patch(SQ_1fishable_0notfishable_for_each_soft_depth_patch); 
     delta_fishable=delta(SQ_1fishable_0notfishable_for_each_soft_depth_patch); 
     %Now run the fleet model in just those patches
     Z=(x0_PPUE+habitat_area_i_fishable.*theta.*(1+gamma.*distance_to_port_for_each_soft_depth_patch_fishable))./(price.*Bi_legal_fishable);
    lambert_z=-exp(-1./Z)./Z;
    %wapr will break if lambert_z is out of range (i.e., <-exp(-1)); so set
    %those values to the edge of the range = -exp(-1)
    lambert_z(lambert_z<-exp(-1))=-exp(-1);
    [lambert_w,~]=wapr(lambert_z); %use stock Matlab function lambertw once get newer version
    deltai_plus_Ei=1./Z+lambert_w';
    deltai_plus_Ei=max([zeros(length(deltai_plus_Ei),1) deltai_plus_Ei],[],2); %get ride of tiny negative results
    Ei_vector_fleet=deltai_plus_Ei-delta_fishable; %remove natural mortality effect
    Ei_vector_fleet=max(Ei_vector_fleet,0); %only count positive harvest patches
    Ei_vector_fleet_all_patches(SQ_1fishable_0notfishable_for_each_soft_depth_patch)=Ei_vector_fleet; %fill in fleet model reults from FISHABLE patches (Ei stays at 0 for rest)
    out=Ei_vector_fleet_all_patches;

end