% plot_map_results
% plot the most recent results on a map
patchmarkersize=14;

%% virgin
metric_to_plot=SSBivirgin;
figure
scatter(lat_lon_habitat_patches(metric_to_plot==0,1),lat_lon_habitat_patches(metric_to_plot==0,2),patchmarkersize,[0.5 0.5 0.5],'s','filled')%
hold on
scatter(lat_lon_habitat_patches(metric_to_plot>0,1),lat_lon_habitat_patches(metric_to_plot>0,2),patchmarkersize,metric_to_plot(metric_to_plot>0),'s','filled')%
% scatter(lat_lon_habitat_patches(SQ_1fishable_0notfishable_for_each_soft_depth_patch==0,1),lat_lon_habitat_patches(SQ_1fishable_0notfishable_for_each_soft_depth_patch==0,2),patchmarkersize,'k','s','filled')
scatter(lat_lon_leftmapedge(1),lat_lon_leftmapedge(2),0.1,'.','k')
colorbar
axis tight
xlabel('Latitude')
ylabel('Longitude')
title('Halibut virgin SSB [kg/km^2/year]')
set(gcf,'color','white');
plot(msp_domain(:,1),msp_domain(:,2),'k','linewidth',coastlinewidth) %coastline
set(gca,'XTickLabel','')
set(gca,'YTickLabel','')
print(strcat(output_figure_dir,'Hal_SSB_virgin_map'))

%% MSY
metric_to_plot=SSBi_wrt_Fsum_msy;
figure
scatter(lat_lon_habitat_patches(metric_to_plot==0,1),lat_lon_habitat_patches(metric_to_plot==0,2),patchmarkersize,[0.5 0.5 0.5],'s','filled')%
hold on
scatter(lat_lon_habitat_patches(metric_to_plot>0,1),lat_lon_habitat_patches(metric_to_plot>0,2),patchmarkersize,metric_to_plot(metric_to_plot>0),'s','filled')%
scatter(lat_lon_leftmapedge(1),lat_lon_leftmapedge(2),0.1,'.','k')
colorbar
axis tight
xlabel('Latitude')
ylabel('Longitude')
title('Halibut SSB at MSY [kg/km^2/year]')
set(gcf,'color','white');
plot(msp_domain(:,1),msp_domain(:,2),'k','linewidth',coastlinewidth) %coastline
set(gca,'XTickLabel','')
set(gca,'YTickLabel','')
print(strcat(output_figure_dir,'Hal_fishery_MSY_SSB_map'))

metric_to_plot=Yi_Fsum_msy;
figure
scatter(lat_lon_habitat_patches(metric_to_plot==0,1),lat_lon_habitat_patches(metric_to_plot==0,2),patchmarkersize,[0.5 0.5 0.5],'s','filled')%
hold on
scatter(lat_lon_habitat_patches(metric_to_plot>0,1),lat_lon_habitat_patches(metric_to_plot>0,2),patchmarkersize,metric_to_plot(metric_to_plot>0),'s','filled')%
scatter(lat_lon_leftmapedge(1),lat_lon_leftmapedge(2),0.1,'.','k')
colorbar
axis tight
xlabel('Latitude')
ylabel('Longitude')
title('Halibut fishery MSY [kg/km^2/year]')
set(gcf,'color','white');
plot(msp_domain(:,1),msp_domain(:,2),'k','linewidth',coastlinewidth) %coastline
set(gca,'XTickLabel','')
set(gca,'YTickLabel','')
print(strcat(output_figure_dir,'Hal_fishery_MSY_yield_map'))

metric_to_plot=Payoffi_Fsum_msy;
figure
scatter(lat_lon_habitat_patches(metric_to_plot==0,1),lat_lon_habitat_patches(metric_to_plot==0,2),patchmarkersize,[0.5 0.5 0.5],'s','filled')%
hold on
scatter(lat_lon_habitat_patches(metric_to_plot>0,1),lat_lon_habitat_patches(metric_to_plot>0,2),patchmarkersize,metric_to_plot(metric_to_plot>0),'s','filled')%
scatter(lat_lon_leftmapedge(1),lat_lon_leftmapedge(2),0.1,'.','k')
colorbar
axis tight
xlabel('Latitude')
ylabel('Longitude')
title('Halibut profit at MSY [$/km^2/year]')
set(gcf,'color','white');
plot(msp_domain(:,1),msp_domain(:,2),'k','linewidth',coastlinewidth) %coastline
set(gca,'XTickLabel','')
set(gca,'YTickLabel','')
print(strcat(output_figure_dir,'Hal_fishery_MSY_profit_map'))
