% plot_CDFWblock_observed_model_data
% plots landings and other data by CDFW block. Plots observed (real) and
% simulated (model) data

%plot observed and simulated results as panels in one figure or as
%individual figures:
plot_individual0_panels1_both2=2; %CHOOSE

shell=NaN(size(unique_block_id));
model_block_sumyield=shell;
% for each CDFW blockid
for index=1:length(unique_block_id)
    patches_in_block=blockid_for_each_soft_depth_patch==unique_block_id(index); %vector with ones=patch in focal block
    model_block_sumyield(index)=sum(Yi(patches_in_block));
end
model_block_sumyield_sum2one=model_block_sumyield./sum(model_block_sumyield); %scale so sums to one

model_sumYi_in_each_block=zeros(size(study_area_polygons_PacCoastFisheryGIS_NUM(:,1)));
for index=1:length(unique_block_id)
    patches_in_block=blockid_for_each_soft_depth_patch==unique_block_id(index); %vector with ones=patch in focal block
    tmp1=study_area_polygons_PacCoastFisheryGIS_NUM(:,block10_id)==unique_block_id(index); %find ALL the patches in the block
    model_sumYi_in_each_block(tmp1)=model_block_sumyield_sum2one(index); %insert summed landings in all the patches in that block
end

%Sum the model Yi values for each CDFW block
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
disp(['gamma=',num2str(gamma),'; psi=',num2str(psi),'; SSE=',num2str(sum(tmp2))])
SSE=sum(tmp2); %goal is to minimize this
% tmp1=mean(block_kglandings_inpatches_sum2one_Rec_Com)-model_sumYi_in_each_block;
% tmp2=tmp1.^2
% SSE_tot=sum(tmp2)
% Rsquare=1-SSE/SSE_tot

%Regress data
[B,BINT,R,RINT,STATS] = regress(block_kglandings_inpatches_sum2one_Rec_Com,[model_sumYi_in_each_block ones(size(model_sumYi_in_each_block))]); disp(['R-square=',num2str(STATS(1))])
% [B,BINT,R,RINT,STATS] = regress(PacCoastFisheryGIS_kglandings_by_block_Com_plus_Rec_sum2one,[model_block_sumyield_sum2one ones(size(model_block_sumyield_sum2one))]); disp(['R-square=',num2str(STATS(1))])
disp('NOTE: R-SQUARE STATISTIC BASED ON LINEAR REGRESSION WITH A NON-ZERO INTERCEPT!!')
%% Individual plots
if plot_individual0_panels1_both2==0 || plot_individual0_panels1_both2==2
    patchmarkersize=14;

    figure
    scatter(lat_lon_msp_domain(:,1),lat_lon_msp_domain(:,2),patchmarkersize,100.*block_kglandings_inpatches_sum2one_Rec_Com,'s','filled') 
    hold on
    scatter(lat_lon_leftmapedge(1),lat_lon_leftmapedge(2),0.1,'.','k')
    colorbar
    axis tight
    xlabel('Latitude')
    ylabel('Longitude')
    title('Halibut average annual total (comm + rec) observed landings [% in SCB]')
    set(gcf,'color','white'); 
    plot(msp_domain(:,1),msp_domain(:,2),'k','linewidth',coastlinewidth) %coastline
%     axis square

    figure
    scatter(lat_lon_msp_domain(:,1),lat_lon_msp_domain(:,2),patchmarkersize,100.*model_sumYi_in_each_block,'s','filled') 
    hold on
    scatter(lat_lon_leftmapedge(1),lat_lon_leftmapedge(2),0.1,'.','k')
    colorbar
    axis tight
    xlabel('Latitude')
    ylabel('Longitude')
    title(['Halibut equil. model landings [% in SCB]; \gamma = ',num2str(gamma),'; \psi = ',num2str(psi),'; SSE = ',num2str(SSE)])
    set(gcf,'color','white'); 
    plot(msp_domain(:,1),msp_domain(:,2),'k','linewidth',coastlinewidth) %coastline
%     axis square

    figure
    plot(100.*block_kglandings_inpatches_sum2one_Rec_Com,100.*model_sumYi_in_each_block,'bo')
    xlabel('Observed landings [% in SCB]')
    ylabel(['Model landings [% in SCB]'])
    title(['\gamma = ',num2str(gamma),'; \psi = ',num2str(psi),'; SSE = ',num2str(SSE),'; R^2 = ',num2str(STATS(1))])
    set(gcf,'color','white'); 
    box off
    hold on
    plot([0 max([100.*block_kglandings_inpatches_sum2one_Rec_Com; 100.*model_sumYi_in_each_block])],[0 max([100.*block_kglandings_inpatches_sum2one_Rec_Com; 100.*model_sumYi_in_each_block])],'k:')
    legend(' Data',' 1:1 line','location','northwest')
end
if plot_individual0_panels1_both2==1 || plot_individual0_panels1_both2==2
    
    patchmarkersize=3;
    
    figure
    subplot(2,2,1)
    scatter(lat_lon_msp_domain(:,1),lat_lon_msp_domain(:,2),patchmarkersize,100.*block_kglandings_inpatches_sum2one_Rec_Com,'s','filled') 
    hold on
    scatter(lat_lon_leftmapedge(1),lat_lon_leftmapedge(2),0.1,'.','k')
    colorbar
    axis tight
%     xlabel('Latitude')
%     ylabel('Longitude')
    title('Observed [CDFW]')
    set(gcf,'color','white'); 
    plot(msp_domain(:,1),msp_domain(:,2),'k','linewidth',coastlinewidth) %coastline
%     axis square
    set(gca,'XTickLabel','')
    set(gca,'YTickLabel','')

    subplot(2,2,2)
    scatter(lat_lon_msp_domain(:,1),lat_lon_msp_domain(:,2),patchmarkersize,100.*model_sumYi_in_each_block,'s','filled') 
    hold on
    scatter(lat_lon_leftmapedge(1),lat_lon_leftmapedge(2),0.1,'.','k')
    colorbar
    axis tight
%     xlabel('Latitude')
%     ylabel('Longitude')
    title(['Simulated [model]'])
    set(gcf,'color','white'); 
    plot(msp_domain(:,1),msp_domain(:,2),'k','linewidth',coastlinewidth) %coastline
%     axis square
    set(gca,'XTickLabel','')
    set(gca,'YTickLabel','')

    suptitle(['Halibut landings [% in SCB]; \gamma = ',num2str(gamma),'; \psi = ',num2str(psi),'; SSE = ',num2str(SSE)])
end

%save model and CDFW data to excel
xlswrite('Model_CDFG_Y_patches',[block_kglandings_inpatches_sum2one_Rec_Com model_sumYi_in_each_block]);
xlswrite('Model_CDFG_Y_blocks',[PacCoastFisheryGIS_kglandings_by_block_Com_plus_Rec_sum2one model_block_sumyield_sum2one]);