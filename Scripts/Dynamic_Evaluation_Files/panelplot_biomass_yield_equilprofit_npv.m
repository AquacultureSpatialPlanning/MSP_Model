% panelplot_biomass_yield_equilprofit_npv
%Generate a panel plot of biomass, yield, equilbrium profit and npv in each patch

patchmarkersize=3;

figure
subplot(2,2,1)
metric_to_plot=sum(Bij,2);
scatter(lat_lon_habitat_patches(metric_to_plot==0,1),lat_lon_habitat_patches(metric_to_plot==0,2),patchmarkersize,[0.5 0.5 0.5],'s','filled')%
hold on
scatter(lat_lon_habitat_patches(metric_to_plot>0,1),lat_lon_habitat_patches(metric_to_plot>0,2),patchmarkersize,metric_to_plot(metric_to_plot>0),'s','filled')%
scatter(lat_lon_leftmapedge(1),lat_lon_leftmapedge(2),0.1,'.','k')
colorbar
axis tight
% xlabel('Latitude')
% ylabel('Longitude')
title('Biomass [kg]')
set(gcf,'color','white'); 
plot(msp_domain(:,1),msp_domain(:,2),'k','linewidth',coastlinewidth) %coastline
set(gca,'XTickLabel','')
set(gca,'YTickLabel','')

subplot(2,2,2)
metric_to_plot=Yi;
scatter(lat_lon_habitat_patches(metric_to_plot==0,1),lat_lon_habitat_patches(metric_to_plot==0,2),patchmarkersize,[0.5 0.5 0.5],'s','filled')%
hold on
scatter(lat_lon_habitat_patches(metric_to_plot>0,1),lat_lon_habitat_patches(metric_to_plot>0,2),patchmarkersize,metric_to_plot(metric_to_plot>0),'s','filled')%
scatter(lat_lon_leftmapedge(1),lat_lon_leftmapedge(2),0.1,'.','k')
colorbar
axis tight
% xlabel('Latitude')
% ylabel('Longitude')
title('Yield [kg]')
set(gcf,'color','white'); 
plot(msp_domain(:,1),msp_domain(:,2),'k','linewidth',coastlinewidth) %coastline
set(gca,'XTickLabel','')
set(gca,'YTickLabel','')

subplot(2,2,3)
metric_to_plot=Payoff;
scatter(lat_lon_habitat_patches(metric_to_plot==0,1),lat_lon_habitat_patches(metric_to_plot==0,2),patchmarkersize,[0.5 0.5 0.5],'s','filled')%
hold on
scatter(lat_lon_habitat_patches(metric_to_plot>0,1),lat_lon_habitat_patches(metric_to_plot>0,2),patchmarkersize,metric_to_plot(metric_to_plot>0),'s','filled')%
scatter(lat_lon_leftmapedge(1),lat_lon_leftmapedge(2),0.1,'.','k')
colorbar
axis tight
% xlabel('Latitude')
% ylabel('Longitude')
title('Equil. profit [$]')
set(gcf,'color','white'); 
plot(msp_domain(:,1),msp_domain(:,2),'k','linewidth',coastlinewidth) %coastline
set(gca,'XTickLabel','')
set(gca,'YTickLabel','')

savefig('Biomass_Yield_Profit_panelplot')
