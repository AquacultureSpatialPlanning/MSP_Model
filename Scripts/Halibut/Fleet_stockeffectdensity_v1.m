function [should_be_zero1]=Fleet_stockeffectdensity_v1(x0_PPUE)
% x0_PPUE = equal ave profit across all patches
% Fleet model with stock effect wrt density

global Bi_legal price delta theta Fsum habitat_area_i gamma distance_to_port_for_each_soft_depth_patch SQ_1fishable_0notfishable_for_each_soft_depth_patch qi

% Ei_vector_fleet=Ei_wrt_PPUE_qi(x0_PPUE,habitat_area_i,theta,price,Bi_legal,delta, gamma, distance_to_port_for_each_soft_depth_patch,SQ_1fishable_0notfishable_for_each_soft_depth_patch,qi); 
Ei_vector_fleet=Ei_wrt_PPUE(x0_PPUE,habitat_area_i,theta,price,Bi_legal,delta, gamma, distance_to_port_for_each_soft_depth_patch,SQ_1fishable_0notfishable_for_each_soft_depth_patch); 

% % %EXPONENTIAL MODEL
% Z=(x0_PPUE+habitat_area_i.*theta)./(price.*Bi_legal);
% lambert_z=-exp(-1./Z)./Z;
% %wapr will break if lambert_z is out of range (i.e., <-exp(-1)); so set
% %those values to the edge of the range = -exp(-1)
% lambert_z(lambert_z<-exp(-1))=-exp(-1);
% [lambert_w nerror]=wapr(lambert_z); %use stock Matlab function lambertw once get newer version
% deltai_plus_Ei=1./Z+lambert_w';
% deltai_plus_Ei=max([zeros(length(deltai_plus_Ei),1) deltai_plus_Ei],[],2); %get ride of tiny negative results
% Ei_vector_fleet=deltai_plus_Ei-delta; %remove natural mortality effect
% Ei_vector_fleet=max(Ei_vector_fleet,0); %only count positive harvest patches

% should_be_zero1=sum(Ei_vector_fleet,1)-Fsum;
should_be_zero1=(sum(Ei_vector_fleet,1)-Fsum)^2;
end
