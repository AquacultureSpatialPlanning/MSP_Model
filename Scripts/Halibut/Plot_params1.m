%% Plot model domain, patch values and demographic parameters

patchmarkersize=14;

%Model domain and area habitat per patch
figure
% scatter(lat_lon_msp_domain(:,1),lat_lon_msp_domain(:,2),patchmarkersize,[0.5 .5 0.5],'s','filled') %all patches=grey filled in squares
hold on
scatter(lat_lon_habitat_patches(:,1),lat_lon_habitat_patches(:,2),patchmarkersize,habitat_area_i,'s','filled')%overlap habitat patches, colored by m2 habitat
scatter(lat_lon_leftmapedge(1),lat_lon_leftmapedge(2),0.1,'.','k')
colorbar
axis tight
xlabel('Latitude')
ylabel('Longitude')
title('Halibut habitat index [m^2]')
set(gcf,'color','white');
plot(msp_domain(:,1),msp_domain(:,2),'k','linewidth',coastlinewidth) %coastline
% axis squarec
print(strcat(output_figure_dir,'Hal_hab_index_map'),'-djpeg')

%Effect of depth on reduced habitat quality
figure
Huagen1990_depth=[10 30 50 70 90];
Haugen1990_freqoccurr=[0.48 0.31 0.11 0.01 0.001];
hold on
plot(Huagen1990_depth,Haugen1990_freqoccurr,'ro')
plot(sort(-depth_i),psi.*log(sort(-depth_i)+1)+1,'b')
plot(Huagen1990_depth,Haugen1990_freqoccurr,'ro')
xlabel('Mean patch depth [m]')
ylabel('Habitat multiplier')
title('Haugen 1990')
set(gcf,'color','white');
box off
ylim([0 1])
legend('Haugen 1990','Model fit')
axis square
print(strcat(output_figure_dir,'Haugen_1990'),'-djpeg')
%Substrate types
patch_color=[1 0 0];
figure
subplot(221)
scatter(lat_lon_msp_domain(hard1_other0_patch_vector,1),lat_lon_msp_domain(hard1_other0_patch_vector,2),patchmarkersize,patch_color,'s','filled')
hold on
scatter(lat_lon_leftmapedge(1),lat_lon_leftmapedge(2),0.1,'.','k')
axis tight
title('Hard')
set(gcf,'color','white');
plot(msp_domain(:,1),msp_domain(:,2),'k','linewidth',coastlinewidth) %coastline
% axis square
subplot(222)
scatter(lat_lon_msp_domain(mixed1_other0_patch_vector,1),lat_lon_msp_domain(mixed1_other0_patch_vector,2),patchmarkersize,patch_color,'s','filled')
hold on
scatter(lat_lon_leftmapedge(1),lat_lon_leftmapedge(2),0.1,'.','k')
axis tight
title('Mixed')
set(gcf,'color','white');
plot(msp_domain(:,1),msp_domain(:,2),'k','linewidth',coastlinewidth) %coastline
% axis square
subplot(223)
scatter(lat_lon_msp_domain(soft1_other0_patch_vector,1),lat_lon_msp_domain(soft1_other0_patch_vector,2),patchmarkersize,patch_color,'s','filled')
hold on
scatter(lat_lon_leftmapedge(1),lat_lon_leftmapedge(2),0.1,'.','k')
axis tight
title('Soft')
set(gcf,'color','white');
plot(msp_domain(:,1),msp_domain(:,2),'k','linewidth',coastlinewidth) %coastline
% axis square
subplot(224)
scatter(lat_lon_msp_domain(unknown1_other0_patch_vector,1),lat_lon_msp_domain(unknown1_other0_patch_vector,2),patchmarkersize,patch_color,'s','filled')
hold on
scatter(lat_lon_leftmapedge(1),lat_lon_leftmapedge(2),0.1,'.','k')
axis tight
title('Unknown')
set(gcf,'color','white');
plot(msp_domain(:,1),msp_domain(:,2),'k','linewidth',coastlinewidth) %coastline
% axis square
print(strcat(output_figure_dir,'Substrate_panel_maps'),'-djpeg')

%Plot of aquaculture developable patches
patch_color=[0 0.5 0];
figure
scatter(lat_lon_msp_domain(dev1_undev0_patch_vector,1),lat_lon_msp_domain(dev1_undev0_patch_vector,2),patchmarkersize,patch_color,'s','filled') %all patches=grey filled in squares
hold on
scatter(lat_lon_leftmapedge(1),lat_lon_leftmapedge(2),0.1,'.','k')
% colorbar
axis tight
xlabel('Latitude')
ylabel('Longitude')
title('Developable for Aquaculture')
set(gcf,'color','white');
plot(msp_domain(:,1),msp_domain(:,2),'k','linewidth',coastlinewidth) %coastline
box on
% axis square
axis tight
print(strcat(output_figure_dir,'Aqua_developable_map'),'-djpeg')

%Plot trawlable patches
patch_color=('red');
figure
scatter(lat_lon_habitat_patches(SQ_1trawling_0notrawling_for_each_soft_detph_patch,1),lat_lon_habitat_patches(SQ_1trawling_0notrawling_for_each_soft_detph_patch,2),patchmarkersize,patch_color,'s','filled') %all patches=grey filled in squares
hold on
scatter(lat_lon_leftmapedge(1),lat_lon_leftmapedge(2),0.1,'.','k')
% axis tight
xlabel('Latitude')
ylabel('Longitude')
title('Trawlable')
set(gcf,'color','white');
plot(msp_domain(:,1),msp_domain(:,2),'k','linewidth',coastlinewidth) %coastline
box on
% axis square
axis tight
print(strcat(output_figure_dir,'Trawlable_map'),'-djpeg')

figure
patch_color=('green');
scatter(lat_lon_habitat_patches(SQ_1fishable_0notfishable_for_each_soft_depth_patch,1),lat_lon_habitat_patches(SQ_1fishable_0notfishable_for_each_soft_depth_patch,2),patchmarkersize,patch_color,'s','filled') %all patches=grey filled in squares
hold on
patch_color=('red');
scatter(lat_lon_habitat_patches(SQ_1fishable_0notfishable_for_each_soft_depth_patch==0,1),lat_lon_habitat_patches(SQ_1fishable_0notfishable_for_each_soft_depth_patch==0,2),patchmarkersize,patch_color,'s','filled')
scatter(lat_lon_leftmapedge(1),lat_lon_leftmapedge(2),0.1,'.','k')
% axis tight
xlabel('Latitude')
ylabel('Longitude')
title('Status quo')
legend('Fishable','Not fishable')
set(gcf,'color','white');
plot(msp_domain(:,1),msp_domain(:,2),'k','linewidth',coastlinewidth) %coastline
box on
% axis square
axis tight
print(strcat(output_figure_dir,'Fishable_notfishable_map'),'-djpeg')

figure
% scatter(lat_lon_habitat_patches(:,1),lat_lon_habitat_patches(:,2),patchmarkersize,patch_color,'s','filled') %all patches=grey filled in squares
hold on
% patch_color=('red');
scatter(lat_lon_habitat_patches(:,1),lat_lon_habitat_patches(:,2),patchmarkersize,depth_i,'s','filled')
scatter(lat_lon_leftmapedge(1),lat_lon_leftmapedge(2),0.1,'.','k')
% axis tight
xlabel('Latitude')
ylabel('Longitude')
title('Depth')
colorbar
set(gcf,'color','white');
plot(msp_domain(:,1),msp_domain(:,2),'k','linewidth',coastlinewidth) %coastline
box on
% axis square
axis tight
print(strcat(output_figure_dir,'Depth_map'),'-djpeg')

%% Plot demographics
figure; plot(0:max_age,Lj_start_at_0age); xlabel('Age [yr]'); ylabel('Length [cm]'); set(gcf,'color',[1 1 1]); box off; grid off; axis square; print(strcat(output_figure_dir,'Age-Length_fcn'),'-djpeg')
figure; plot(0:max_age,Wj_start_at_0age); xlabel('Age [yr]'); ylabel('Weight [kg]'); set(gcf,'color',[1 1 1]); box off; grid off; axis square; print(strcat(output_figure_dir,'Age-Weight fcn'),'-djpeg')
figure; surface(Dii); title('Larval dispersal'); xlabel('Destination'); ylabel('Source'); set(gcf,'color',[1 1 1]); colorbar; shading flat; axis tight; axis square; print(strcat(output_figure_dir,'Larval_dispersal_matrix'),'-djpeg')
% figure; surface(Mii); title('Adult movement'); xlabel('Destination'); ylabel('Source'); set(gcf,'color',[1 1 1]); colorbar; shading flat; axis tight; print(strcat(output_figure_dir,'Adult_movement_matrix'),'-djpeg')

%% Economics
figure
hist(1+gamma.*distance_to_port_for_each_soft_depth_patch)
ylabel('Number of patches')
xlabel('Cost factor increase due to distance from port')
set(gcf,'color',[1 1 1]);
box off
print(strcat(output_figure_dir,'Hist_patches_costfactorincreasewrtportdist'),'-djpeg')

figure
plot(distance_to_port_for_each_soft_depth_patch./1000, (1+gamma.*distance_to_port_for_each_soft_depth_patch),'bo')
xlabel('Patch distance from port [km]')
ylabel('Cost factor increase')
set(gcf,'color',[1 1 1]);
box off
print(strcat(output_figure_dir,'Costfactorincrease_wrt_patchportdist'),'-djpeg')

figure
hold on
scatter(lat_lon_habitat_patches(:,1),lat_lon_habitat_patches(:,2),patchmarkersize,distance_to_port_for_each_soft_depth_patch./1000,'s','filled')
scatter(lat_lon_leftmapedge(1),lat_lon_leftmapedge(2),0.1,'.','k')
% axis tight
xlabel('Latitude')
ylabel('Longitude')
title('Distance from port [km]')
colorbar
set(gcf,'color','white');
plot(msp_domain(:,1),msp_domain(:,2),'k','linewidth',coastlinewidth) %coastline
box on
axis tight
print(strcat(output_figure_dir,'Map_portdist'),'-djpeg')

figure
hold on
scatter(lat_lon_habitat_patches(:,1),lat_lon_habitat_patches(:,2),patchmarkersize,(1+gamma.*distance_to_port_for_each_soft_depth_patch),'s','filled')
scatter(lat_lon_leftmapedge(1),lat_lon_leftmapedge(2),0.1,'.','k')
% axis tight
xlabel('Latitude')
ylabel('Longitude')
title('Cost factor increase due to distance from port')
colorbar
set(gcf,'color','white');
plot(msp_domain(:,1),msp_domain(:,2),'k','linewidth',coastlinewidth) %coastline
box on
axis tight
print(strcat(output_figure_dir,'Map_costfactorincrease'),'-djpeg')
