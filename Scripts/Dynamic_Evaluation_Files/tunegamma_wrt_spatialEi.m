function [out]=tunegamma_wrt_spatialEi(gamma_input)

global habitat_area_i numpatches max_age Bij age_mature Dii K age_legal Ei delta Nij Mii age_move Wij TargetBiomass beta_i...
    Rmax CRgoal Nij_initial theta price Fsum_msy_guess SSB_msy_target Theta_Density_prop Fsum alphaCR gamma distance_to_port_for_each_soft_depth_patch unique_block_id...
    blockid_for_each_soft_depth_patch PacCoastFisheryGIS_kglandings_by_block_Com_plus_Rec_sum2one yield_scalar SQ_1fishable_0notfishable_for_each_soft_depth_patch qi omega nearestdist_to_outfall...
    tol_simulation T x0_PPUE_initial discount_rate

gamma=gamma_input; %consider gamma distance from port effect
BasicSpatialHarvestModel_COREv1_fleet_TS %run model wrt the gamma

shell=NaN(size(unique_block_id));
model_block_sumyield=shell;
% for each CDFW blockid
for index=1:length(unique_block_id)
    patches_in_block=blockid_for_each_soft_depth_patch==unique_block_id(index); %vector with ones=patch in focal block
    model_block_sumyield(index)=sum(Yi(patches_in_block));
end
model_block_sumyield_sum2one=model_block_sumyield./sum(model_block_sumyield); %scale so sums to one
%calculate the sum of the squared difference between the model and CDFW block landings
tmp1=PacCoastFisheryGIS_kglandings_by_block_Com_plus_Rec_sum2one-model_block_sumyield_sum2one; %diff between model and CDFW block landings data
tmp2=tmp1.^2; %diff squared
disp(['gamma_input=',num2str(gamma_input),'; SSE=',num2str(sum(tmp2))])
out=sum(tmp2); %goal is to minimize this
end

% %% Evaluate effect of gamma manually - comment out normally
% % close all
% % gamma_range=linspace(0,0.005,10);
% % gamma_range=linspace(0,0.01,10);
% % gamma_range=[0:0.001:0.01];
% % gamma_range=[0:0.0001:0.001];
% gamma_range=linspace(0,2*2.7379e-05,10);
% SSE_wrt_gamma=NaN(size(gamma_range));
% gamma_index=0;
% for gamma=gamma_range
%     gamma_index=gamma_index+1;
%     disp(['gamma=',num2str(gamma)])
%     
%     BasicSpatialHarvestModel_COREv1_fleet %run model wrt the gamma
% 
%     shell=NaN(size(unique_block_id));
%     model_block_sumyield=shell;
%     % for each CDFW blockid
%     for index=1:length(unique_block_id)
%         patches_in_block=blockid_for_each_soft_depth_patch==unique_block_id(index); %vector with ones=patch in focal block
%         model_block_sumyield(index)=sum(Yi(patches_in_block));
%     end
%     model_block_sumyield_sum2one=model_block_sumyield./sum(model_block_sumyield); %scale so sums to one
%     %calculate the sum of the squared difference between the model and CDFW block landings
%     tmp1=PacCoastFisheryGIS_kglandings_by_block_Com_plus_Rec_sum2one-model_block_sumyield_sum2one; %diff between model and CDFW block landings data
%     tmp2=tmp1.^2; %diff squared
%     disp(['gamma=',num2str(gamma),'; SSE=',num2str(sum(tmp2))])
%     out=sum(tmp2); %goal is to minimize this
%     SSE_wrt_gamma(gamma_index)=out;
% end
% figure
% plot(gamma_range,SSE_wrt_gamma)
% xlabel('\gamma')
% ylabel('SSE') 
% set(gcf,'color','white'); 
% box off
% 
%% plot spatial distribution of landings - comment out normally
% %Calculate SSE
% BasicSpatialHarvestModel_COREv1_fleet %run model wrt the gamma
% 
% shell=NaN(size(unique_block_id));
% model_block_sumyield=shell;
% % for each CDFW blockid
% for index=1:length(unique_block_id)
%     patches_in_block=blockid_for_each_soft_depth_patch==unique_block_id(index); %vector with ones=patch in focal block
%     model_block_sumyield(index)=sum(Yi(patches_in_block));
% end
% model_block_sumyield_sum2one=model_block_sumyield./sum(model_block_sumyield); %scale so sums to one
% 
% model_sumYi_in_each_block=zeros(size(study_area_polygons_PacCoastFisheryGIS_NUM(:,1)));
% for index=1:length(unique_block_id)
%     patches_in_block=blockid_for_each_soft_depth_patch==unique_block_id(index); %vector with ones=patch in focal block
%     tmp1=study_area_polygons_PacCoastFisheryGIS_NUM(:,block10_id)==unique_block_id(index); %find ALL the patches in the block
%     model_sumYi_in_each_block(tmp1)=model_block_sumyield_sum2one(index); %insert summed landings in all the patches in that block
% end
% 
% %Sum the model Yi values for each CDFW block
% shell=NaN(size(unique_block_id));
% model_block_sumyield=shell;
% % for each CDFW blockid
% for index=1:length(unique_block_id)
%     patches_in_block=blockid_for_each_soft_depth_patch==unique_block_id(index); %vector with ones=patch in focal block
%     model_block_sumyield(index)=sum(Yi(patches_in_block));
% end
% model_block_sumyield_sum2one=model_block_sumyield./sum(model_block_sumyield); %scale so sums to one
% %calculate the sum of the squared difference between the model and CDFW block landings
% tmp1=PacCoastFisheryGIS_kglandings_by_block_Com_plus_Rec_sum2one-model_block_sumyield_sum2one; %diff between model and CDFW block landings data
% tmp2=tmp1.^2; %diff squared
% disp(['gamma=',num2str(gamma),'; SSE=',num2str(sum(tmp2))])
% SSE=sum(tmp2); %goal is to minimize this
% 
% figure
% scatter(lat_lon_msp_domain(:,1),lat_lon_msp_domain(:,2),patchmarkersize,100.*block_kglandings_inpatches_sum2one_Rec_Com,'s','filled') 
% hold on
% scatter(lat_lon_leftmapedge(1),lat_lon_leftmapedge(2),0.1,'.','k')
% colorbar
% axis tight
% xlabel('Latitude')
% ylabel('Longitude')
% title('Halibut average annual total (comm + rec) landings [% in SCB]')
% set(gcf,'color','white'); 
% plot(msp_domain(:,1),msp_domain(:,2),'k','linewidth',coastlinewidth) %coastline
% axis square
% 
% figure
% scatter(lat_lon_msp_domain(:,1),lat_lon_msp_domain(:,2),patchmarkersize,100.*model_sumYi_in_each_block,'s','filled') 
% hold on
% scatter(lat_lon_leftmapedge(1),lat_lon_leftmapedge(2),0.1,'.','k')
% colorbar
% axis tight
% xlabel('Latitude')
% ylabel('Longitude')
% title(['Halibut equilibrium model landings [% in SCB]; \gamma = ',num2str(gamma),'; SSE = ',num2str(SSE)])
% set(gcf,'color','white'); 
% plot(msp_domain(:,1),msp_domain(:,2),'k','linewidth',coastlinewidth) %coastline
% axis square
% 
% %Regress data
% [B,BINT,R,RINT,STATS] = regress(block_kglandings_inpatches_sum2one_Rec_Com,model_sumYi_in_each_block);
% 
% figure
% plot(100.*block_kglandings_inpatches_sum2one_Rec_Com,100.*model_sumYi_in_each_block,'bo')
% xlabel('Comm + rec landings [% in SCB]')
% ylabel(['Model landings [% in SCB]'])
% title(['\gamma = ',num2str(gamma),'; SSE = ',num2str(SSE),' R^2 = ',num2str(STATS(1))])
% set(gcf,'color','white'); 
% box off
% hold on
% plot([0 max([100.*block_kglandings_inpatches_sum2one_Rec_Com; 100.*model_sumYi_in_each_block])],[0 max([100.*block_kglandings_inpatches_sum2one_Rec_Com; 100.*model_sumYi_in_each_block])],'k:')
% 

