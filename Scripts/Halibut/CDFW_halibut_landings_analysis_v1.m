%Load and analyze PacCoastFisheryGIS data from CDFW
%Email communications:
% Serpa, Paulo@Wildlife <Paulo.Serpa@wildlife.ca.gov>
% Paulo Serpa
% Geographic Information Systems
% Statewide MPA Management Project
% California Department of Fish and Wildlife | Marine Region
% 20 Lower Ragsdale Drive #100
% Monterey, CA 93940
% (831) 649 - 7143
% http://www.dfg.ca.gov/marine/gis/
% Also on the email:
%  "Aseltine-Neilson, Debbie@Wildlife" <Debbie.Aseltine-Neilson@wildlife.ca.gov>,
%  "Eres, Joann@Wildlife" <Joann.Eres@wildlife.ca.gov> 

clear all
close all

set(0,'DefaultTextFontSize',35)
set(0,'DefaultAxesFontSize',35)
set(0,'DefaultLineLineWidth',6)
set(0,'DefaultAxesLineWidth',2)
set(0,'DefaultSurfaceLineWidth',2)

years_data_to_collect=2004:2100; %just the most recent years corresponding with Maunder et al. 2011 and when the stock and yield was steady
%% California halibut (Paralichthys californicus) ---------
%Load GIS data from Carrie Kappel
% [study_area_polygons_PacCoastFisheryGIS_NUM,study_area_polygons_PacCoastFisheryGIS_TXT,study_area_polygons_PacCoastFisheryGIS_RAW]=xlsread('SeaGrant_allcells_FinalJuly2014.xlsx');
[study_area_polygons_PacCoastFisheryGIS_NUM,study_area_polygons_PacCoastFisheryGIS_TXT,study_area_polygons_PacCoastFisheryGIS_RAW]=xlsread('SeaGrant_allFINALdata_complete_2015.xlsx');
col_lat=2; %x-axis = latitude
col_lon=3; %y-axis = longitude
block_id_col=15; %row with CDFW blocks corresponding with our study patches

%ID the unique block IDs, and get rid of NaNs
unique_block_id=unique(study_area_polygons_PacCoastFisheryGIS_NUM(:,block_id_col));
unique_block_id(isnan(unique_block_id))=[];
disp(['Our 1-km2 patches in the SCB lie within ',num2str(length(unique_block_id)),' unique CDFW blocks'])

%% Load PacCoastFisheryGIS data: Commercial
[Commercial_Hal_NUM,Commercial_Hal_TXT,Commercial_Hal_RAW]=xlsread('Commercial_Hal.xlsx');
year_col_Com=2;
month_col_Com=3;
block_code_col_Com=4;
pounds_col_Com=9;

%Filter so just looking at data within years_data_to_collect
Commercial_Hal_NUM=Commercial_Hal_NUM(ismember(Commercial_Hal_NUM(:,year_col_Com),years_data_to_collect),:);

%It doeesn't look like the block_codes are integers, so pull them and round
%them so can match them with the integers in unique_block_id
block_code_Com=round(Commercial_Hal_NUM(:,block_code_col_Com));

max_years_possible=length(years_data_to_collect); %50;
lbslanding_by_block_year_Com=NaN(length(unique_block_id),max_years_possible);

%for each block_id in 'unique_block_id'...
for indexb=1:length(unique_block_id)
    block_id=unique_block_id(indexb);
%      find the corresponding rows of block codes in 'Commercial_Hal_NUM'
    Commercial_Hal_NUM_blockrows=Commercial_Hal_NUM(block_code_Com==block_id,:);
    %for each year
    unique_year=unique(Commercial_Hal_NUM_blockrows(:,year_col_Com));
    for indexy=1:length(unique_year)
        year=unique_year(indexy);
        %sum all the months of landings
        Commercial_Hal_NUM_blockrows_year=Commercial_Hal_NUM_blockrows(Commercial_Hal_NUM_blockrows(:,year_col_Com)==year,:);
        lbslanding_by_block_year_Com(indexb,indexy)=nansum(Commercial_Hal_NUM_blockrows_year(:,pounds_col_Com));
    end
end
    
%average landing in each block across years
lbslanding_by_block_Com=nanmean(lbslanding_by_block_year_Com,2);
%if lbslanding_by_block_Com=[] then no landings in that block for the years of the analysis, so set to zero
lbslanding_by_block_Com(isnan(lbslanding_by_block_Com))=0;
%convert units
lb_per_kg=2.20462;
PacCoastFisheryGIS_kglandings_by_block_Com=lbslanding_by_block_Com./lb_per_kg;
%calculate proportional values that sum to 1
PacCoastFisheryGIS_kglandings_by_block_sum2one_Com=PacCoastFisheryGIS_kglandings_by_block_Com./sum(PacCoastFisheryGIS_kglandings_by_block_Com);

% assign landings block data to patches in each block
shell=NaN(length(study_area_polygons_PacCoastFisheryGIS_NUM(:,1)),1);
block_kglandings_inpatches_Com=shell;
block_kglandings_inpatches_sum2one_Com=shell;
%for each block..
for indexb=1:length(unique_block_id)
    block_id=unique_block_id(indexb);
    
    %find what patches are in that block, and assign the block landings
    %data to those patch rows
    block_kglandings_inpatches_Com(study_area_polygons_PacCoastFisheryGIS_NUM(:,block_id_col)==block_id)=PacCoastFisheryGIS_kglandings_by_block_Com(indexb);
    block_kglandings_inpatches_sum2one_Com(study_area_polygons_PacCoastFisheryGIS_NUM(:,block_id_col)==block_id)=PacCoastFisheryGIS_kglandings_by_block_sum2one_Com(indexb);
end
    
% plot the patches wrt the block landings data
%coastline
msp_domain = textread('msp_domain.dat');
map_border_buffer_prop=0.001; %proportion shoreline map larger than edge patches
%Because the SCB coast does not cross the left border of the map, the
%patches end up being the left-most point (not good for display). So,
%create a point that is slightly left of the patches and have that be the edge :)
lat_lon_msp_domain=study_area_polygons_PacCoastFisheryGIS_NUM(:,[col_lat col_lon]); %matrix of lat and lon coordinates to every patch in the MSP study area
lat_lon_leftmapedge=[min(lat_lon_msp_domain(:,1))*(1+map_border_buffer_prop) max(lat_lon_msp_domain(:,2))*(1+map_border_buffer_prop)];
%set mapping symbols and line dimensions (NEED TO WORK ON THESE NUMBERS!)
coastlinewidth=2; %thickness of coastline line on maps
patchmarkersize=9; %size of patch squares on maps (set this so that the patches to "fit" together like a perfect grid)

figure
scatter(lat_lon_msp_domain(:,1),lat_lon_msp_domain(:,2),patchmarkersize,block_kglandings_inpatches_Com,'s','filled') %all patches=grey filled in squares
% hold on
% scatter(lat_lon_habitat_patches(:,1),lat_lon_habitat_patches(:,2),patchmarkersize,block_kglandings_inpatches,'s','filled')%overlap habitat patches, colored by m2 habitat
hold on
scatter(lat_lon_leftmapedge(1),lat_lon_leftmapedge(2),0.1,'.','k')
colorbar
axis tight
xlabel('Latitude')
ylabel('Longitude')
title('Halibut average annual commercial landings [kg]')
set(gcf,'color','white'); 
plot(msp_domain(:,1),msp_domain(:,2),'k','linewidth',coastlinewidth) %coastline
axis square
savefig('Hal_comm_kg_landings')

figure
scatter(lat_lon_msp_domain(:,1),lat_lon_msp_domain(:,2),patchmarkersize,100.*block_kglandings_inpatches_sum2one_Com,'s','filled') %all patches=grey filled in squares
% hold on
% scatter(lat_lon_habitat_patches(:,1),lat_lon_habitat_patches(:,2),patchmarkersize,block_kglandings_inpatches,'s','filled')%overlap habitat patches, colored by m2 habitat
hold on
scatter(lat_lon_leftmapedge(1),lat_lon_leftmapedge(2),0.1,'.','k')
colorbar
axis tight
xlabel('Latitude')
ylabel('Longitude')
title('Halibut average annual commercial landings [% in SCB]')
set(gcf,'color','white'); 
plot(msp_domain(:,1),msp_domain(:,2),'k','linewidth',coastlinewidth) %coastline
axis square
savefig('Hal_comm_percent_landings')

%% Load PacCoastFisheryGIS data: Recreational landings data
[Recreational_Hal_NUM,Recreational_Hal_TXT,Recreational_Hal_RAW]=xlsread('Recreational_Hal.xlsx');
SumOf_Kept_Rec_col=3;
Block_Rec_col=4;
Year_Rec_col=5;
% Month_Rec_col=6;

%Filter so just looking at data within years_data_to_collect
Recreational_Hal_NUM=Recreational_Hal_NUM(ismember(Recreational_Hal_NUM(:,Year_Rec_col),years_data_to_collect),:);

%It doeesn't look like the block_codes are integers, so pull them and round
%them so can match them with the integers in unique_block_id
block_code_Rec=round(Recreational_Hal_NUM(:,Block_Rec_col));

max_years_possible=length(years_data_to_collect); %50;
SumOfKept_by_block_year=NaN(length(unique_block_id),max_years_possible);

%for each block_id in 'unique_block_id'...
for indexb=1:length(unique_block_id)
    block_id=unique_block_id(indexb);
%      find the corresponding rows of block codes in 'Commercial_Hal_NUM'
    Recreational_Hal_NUM_blockrows=Recreational_Hal_NUM(block_code_Rec==block_id,:);
    %for each year
    unique_year=unique(Recreational_Hal_NUM_blockrows(:,Year_Rec_col));
    for indexy=1:length(unique_year)
        year=unique_year(indexy);
        %sum all the months of landings
        Recreational_Hal_NUM_blockrows_year=Recreational_Hal_NUM_blockrows(Recreational_Hal_NUM_blockrows(:,Year_Rec_col)==year,:);
        SumOfKept_by_block_year(indexb,indexy)=nansum(Recreational_Hal_NUM_blockrows_year(:,SumOf_Kept_Rec_col));
    end
end
    
%average landing in each block across years
SumOfKept_by_block=nanmean(SumOfKept_by_block_year,2);
SumOfKept_by_block(isnan(SumOfKept_by_block))=0; %if no fish landed in that block, then set NaN to zero
kg_per_keptfish=2.7686; %3; %the approx median of the weight frequency distribution of 
% recreationally (CPFV and other) landed and retained halibut (Maunder et al. 2011 CA halibut stock assessment Fig. B2.8.2)
PacCoastFisheryGIS_kglandings_by_block_Rec=SumOfKept_by_block.*kg_per_keptfish;
PacCoastFisheryGIS_kglandings_by_block_sum2one_Rec=PacCoastFisheryGIS_kglandings_by_block_Rec./sum(PacCoastFisheryGIS_kglandings_by_block_Rec);

% assign landings block data to patches in each block
shell=NaN(length(study_area_polygons_PacCoastFisheryGIS_NUM(:,1)),1);
block_kglandings_inpatches_Rec=shell;
block_kglandings_inpatches_sum2one_Rec=shell;
%for each block..
for indexb=1:length(unique_block_id)
    block_id=unique_block_id(indexb);
    
    %find what patches are in that block, and assign the block landings
    %data to those patch rows
    block_kglandings_inpatches_Rec(study_area_polygons_PacCoastFisheryGIS_NUM(:,block_id_col)==block_id)=PacCoastFisheryGIS_kglandings_by_block_Rec(indexb);
    block_kglandings_inpatches_sum2one_Rec(study_area_polygons_PacCoastFisheryGIS_NUM(:,block_id_col)==block_id)=PacCoastFisheryGIS_kglandings_by_block_sum2one_Rec(indexb);
end

% plot the patches wrt the block landings data
figure
scatter(lat_lon_msp_domain(:,1),lat_lon_msp_domain(:,2),patchmarkersize,block_kglandings_inpatches_Rec,'s','filled') %all patches=grey filled in squares
% hold on
% scatter(lat_lon_habitat_patches(:,1),lat_lon_habitat_patches(:,2),patchmarkersize,block_kglandings_inpatches,'s','filled')%overlap habitat patches, colored by m2 habitat
hold on
scatter(lat_lon_leftmapedge(1),lat_lon_leftmapedge(2),0.1,'.','k')
colorbar
axis tight
xlabel('Latitude')
ylabel('Longitude')
title('Halibut average annual recreational landings [kg]')
set(gcf,'color','white'); 
plot(msp_domain(:,1),msp_domain(:,2),'k','linewidth',coastlinewidth) %coastline
axis square
savefig('Hal_rec_kg_landings')

figure
scatter(lat_lon_msp_domain(:,1),lat_lon_msp_domain(:,2),patchmarkersize,100.*block_kglandings_inpatches_sum2one_Rec,'s','filled') %all patches=grey filled in squares
% hold on
% scatter(lat_lon_habitat_patches(:,1),lat_lon_habitat_patches(:,2),patchmarkersize,block_kglandings_inpatches,'s','filled')%overlap habitat patches, colored by m2 habitat
hold on
scatter(lat_lon_leftmapedge(1),lat_lon_leftmapedge(2),0.1,'.','k')
colorbar
axis tight
xlabel('Latitude')
ylabel('Longitude')
title('Halibut average annual recreational landings [% in SCB]')
set(gcf,'color','white'); 
plot(msp_domain(:,1),msp_domain(:,2),'k','linewidth',coastlinewidth) %coastline
axis square
savefig('Hal_rec_percent_landings')

%% Combine commercial and recreational landings data
% disp(['Commercial landings is ~',num2str(sum(PacCoastFisheryGIS_kglandings_by_block_Com)/sum(PacCoastFisheryGIS_kglandings_by_block_Rec)),' times greater than recrational landings'])
disp(['Commercial landings is ',num2str(sum(PacCoastFisheryGIS_kglandings_by_block_Com)/sum(PacCoastFisheryGIS_kglandings_by_block_Rec+PacCoastFisheryGIS_kglandings_by_block_Com)),' proportion of total landings'])
%Adjust landings so that match Maunder et al. 2011 Comm:Rec proportion??
adjust_landings=1;
if adjust_landings==1
    Com_prop_total_Maunderetal2011=0.636905092; %0.615601284; %comm:(comm+rec) landings from Maunder et al 2011 (see Table B1.6.1 in Maunder et al. 2011 CA halibut  stock asessment figures.xlsx)
    Rec_multiplier=(sum(PacCoastFisheryGIS_kglandings_by_block_Com)*(1/Com_prop_total_Maunderetal2011)-sum(PacCoastFisheryGIS_kglandings_by_block_Com))/sum(PacCoastFisheryGIS_kglandings_by_block_Rec);
    PacCoastFisheryGIS_kglandings_by_block_Rec=PacCoastFisheryGIS_kglandings_by_block_Rec.*Rec_multiplier;
    disp(['ADJUSTED:Commercial landings is ',num2str(sum(PacCoastFisheryGIS_kglandings_by_block_Com)/sum(PacCoastFisheryGIS_kglandings_by_block_Rec+PacCoastFisheryGIS_kglandings_by_block_Com)),' proportion of total landings'])
end

PacCoastFisheryGIS_kglandings_by_block_Com_plus_Rec=PacCoastFisheryGIS_kglandings_by_block_Com+PacCoastFisheryGIS_kglandings_by_block_Rec;
PacCoastFisheryGIS_kglandings_by_block_Com_plus_Rec_sum2one=PacCoastFisheryGIS_kglandings_by_block_Com_plus_Rec./sum(PacCoastFisheryGIS_kglandings_by_block_Com_plus_Rec);

% assign landings block data to patches in each block
shell=NaN(length(study_area_polygons_PacCoastFisheryGIS_NUM(:,1)),1);
block_kglandings_inpatches_Rec_Com=shell;
block_kglandings_inpatches_sum2one_Rec_Com=shell;
%for each block..
for indexb=1:length(unique_block_id)
    block_id=unique_block_id(indexb);
    
    %find what patches are in that block, and assign the block landings
    %data to those patch rows
    block_kglandings_inpatches_Rec_Com(study_area_polygons_PacCoastFisheryGIS_NUM(:,block_id_col)==block_id)=PacCoastFisheryGIS_kglandings_by_block_Com_plus_Rec(indexb);
    block_kglandings_inpatches_sum2one_Rec_Com(study_area_polygons_PacCoastFisheryGIS_NUM(:,block_id_col)==block_id)=PacCoastFisheryGIS_kglandings_by_block_Com_plus_Rec_sum2one(indexb);
end

% plot the patches wrt the block landings data
figure
scatter(lat_lon_msp_domain(:,1),lat_lon_msp_domain(:,2),patchmarkersize,block_kglandings_inpatches_Rec_Com,'s','filled') %all patches=grey filled in squares
% hold on
% scatter(lat_lon_habitat_patches(:,1),lat_lon_habitat_patches(:,2),patchmarkersize,block_kglandings_inpatches,'s','filled')%overlap habitat patches, colored by m2 habitat
hold on
scatter(lat_lon_leftmapedge(1),lat_lon_leftmapedge(2),0.1,'.','k')
colorbar
axis tight
xlabel('Latitude')
ylabel('Longitude')
title('Halibut average annual total (comm + rec) landings [kg]')
set(gcf,'color','white'); 
plot(msp_domain(:,1),msp_domain(:,2),'k','linewidth',coastlinewidth) %coastline
axis square
savefig('Hal_total_kg_landings')

figure
scatter(lat_lon_msp_domain(:,1),lat_lon_msp_domain(:,2),patchmarkersize,100.*block_kglandings_inpatches_sum2one_Rec_Com,'s','filled') %all patches=grey filled in squares
% hold on
% scatter(lat_lon_habitat_patches(:,1),lat_lon_habitat_patches(:,2),patchmarkersize,block_kglandings_inpatches,'s','filled')%overlap habitat patches, colored by m2 habitat
hold on
scatter(lat_lon_leftmapedge(1),lat_lon_leftmapedge(2),0.1,'.','k')
colorbar
axis tight
xlabel('Latitude')
ylabel('Longitude')
title('Halibut average annual total (comm + rec) landings [% in SCB]')
set(gcf,'color','white'); 
plot(msp_domain(:,1),msp_domain(:,2),'k','linewidth',coastlinewidth) %coastline
axis square
savefig('Hal_total_percent_landings')

%% Save data needed for the analysis
save('PacCoastFisheryGIS_kglandings_by_block_Com_plus_Rec_sum2one','PacCoastFisheryGIS_kglandings_by_block_Com_plus_Rec_sum2one','unique_block_id','block_kglandings_inpatches_sum2one_Rec_Com','PacCoastFisheryGIS_kglandings_by_block_Com_plus_Rec')
