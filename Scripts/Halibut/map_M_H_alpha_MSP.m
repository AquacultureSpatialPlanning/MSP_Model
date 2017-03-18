metric_to_plot=Yi;
figure
scatter(lat_lon_habitat_patches(metric_to_plot==0,1),lat_lon_habitat_patches(metric_to_plot==0,2),patchmarkersize,[0.5 0.5 0.5],'s','filled')%
hold on
scatter(lat_lon_habitat_patches(metric_to_plot>0,1),lat_lon_habitat_patches(metric_to_plot>0,2),patchmarkersize,metric_to_plot(metric_to_plot>0),'s','filled')%
scatter(lat_lon_positivevaluemussel(X==1,1),lat_lon_positivevaluemussel(X==1,2),patchmarkersize,[0 0 0],'s','filled')%
scatter(lat_lon_leftmapedge(1),lat_lon_leftmapedge(2),0.1,'.','k')
colorbar
axis tight
xlabel('Latitude')
ylabel('Longitude')
title(['Yield (color) and aquaculture (black), given \alpha = ',num2str(alpha)])
set(gcf,'color','white'); 
plot(msp_domain(:,1),msp_domain(:,2),'k','linewidth',coastlinewidth) %coastline
set(gca,'XTickLabel','')
set(gca,'YTickLabel','')


