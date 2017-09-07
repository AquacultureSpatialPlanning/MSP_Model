%% Set model parameters
%Nomenclature of subscripts to parameters used below
% i = patch = rows
% j = age = column

%% California halibut (Paralichthys californicu) parameter values ---------
%Indicate what columns show what data
col_target_fid=1;
col_lat=2; %x-axis = latitude
col_lon=3; %y-axis = longitude
col_class=4; %undev=undevelopable for aquaculture; dev=develaple for aquaculture
 %col # indicating area of habitat in each cell
col_hard_hab_area=5;
col_mixed_hab_area=6;
col_soft_hab_area=7;
col_unknown_hab_area=8;
col_mean_depth=12;
col_distance_to_port=14; %meters to the nearest port
block10_id=15; %CDFW fishing block
col_MPA=19;
col_PLATFORM=20;
col_NAVIGATION=21;
col_ANCHORAGE=22;
col_MILITARY=23;
col_OTHER=24;
col_Trawlable=27; %Whether or not patch is trawlable. Yes if in CA Halibut Trawling Ground (CHTG) or in Federal waters (i.e., outside State waters)
col_rec_fish=28; %recreational fishing allowed
col_mussel_value=32; %mussel values (those in parentheses are negative)
%NOTE: the trawling spatial regulations were implemented as of August 12, 2008
% Carrie confirmed with Paulo and the final rule making that the status is actually correct on those subareas.
% Turns out that though all four were proposed for closure, the fish and game commission opted to only close one of them - area B.
% Here's  a link to the final regulatory action, which postdates the report that you shared: http://www.fgc.ca.gov/regulations/2008/124fsor.pdf.
% You can see the whole policy progression here: http://www.fgc.ca.gov/regulations/2008/#124.
col_outfall=29; %Cells that intersect a wastewater outfall

%Load GIS data from Carrie Kappel
if runsetparams_1full==1
%     [study_area_polygons_PacCoastFisheryGIS_NUM,study_area_polygons_PacCoastFisheryGIS_TXT,study_area_polygons_PacCoastFisheryGIS_RAW]=xlsread('SeaGrant_allcells_FinalJuly2014.xlsx');
[study_area_polygons_PacCoastFisheryGIS_NUM,study_area_polygons_PacCoastFisheryGIS_TXT,study_area_polygons_PacCoastFisheryGIS_RAW]=xlsread(strcat(input_data_dir,'SeaGrant_allFINALdata_complete_2015.xlsx'));
    study_area_polygons_PacCoastFisheryGIS_TXT_noheader=study_area_polygons_PacCoastFisheryGIS_TXT(2:end,:); %remove header from text matrix
end

target_fid_fulldomain=study_area_polygons_PacCoastFisheryGIS_NUM(:,col_target_fid);
lat_lon_msp_domain=study_area_polygons_PacCoastFisheryGIS_NUM(:,[col_lat col_lon]); %matrix of lat and lon coordinates to every patch in the MSP study area

%%

%create matix of values for just halibut habitat (soft and mixed)) cells
%filter rows wrt positive habitat  within depth range:
min_percent_soft_mixed=0;
hard=study_area_polygons_PacCoastFisheryGIS_NUM(:,col_hard_hab_area);
hard(isnan(hard))=0; %isnans represent zero of this habitat, so set to zero
soft=study_area_polygons_PacCoastFisheryGIS_NUM(:,col_soft_hab_area);
soft(isnan(soft))=0; %isnans represent zero of this habitat, so set to zero
mixed=study_area_polygons_PacCoastFisheryGIS_NUM(:,col_mixed_hab_area);
mixed(isnan(mixed))=0;%isnans represent zero of this habitat, so set to zero

%If zero CFG landings in a block, then set habitat of patches in that block to zero
%Empirical CDFW spatial block landings data
load(strcat(input_data_dir,'PacCoastFisheryGIS_kglandings_by_block_Com_plus_Rec_sum2one'),'block_kglandings_inpatches_sum2one_Rec_Com')
soft(block_kglandings_inpatches_sum2one_Rec_Com==0)=0;
mixed(block_kglandings_inpatches_sum2one_Rec_Com==0)=0;

%Set all habitat north of 34.6 latitude equal to zero
max_lat_hab=34.6;
soft(lat_lon_msp_domain(:,2)>max_lat_hab)=0;
mixed(lat_lon_msp_domain(:,2)>max_lat_hab)=0;

unknown=study_area_polygons_PacCoastFisheryGIS_NUM(:,col_unknown_hab_area);
unknown(isnan(unknown))=0;%isnans represent zero of this habitat, so set to zero
total_habitat=hard+soft+mixed+unknown; %note does not always sum to 1, so need to calculate total for each patch
percent_soft_mixed=(soft+mixed)./total_habitat;
filter_habitat_rows=percent_soft_mixed>min_percent_soft_mixed&study_area_polygons_PacCoastFisheryGIS_NUM(:,col_mean_depth)>halibut_max_depth; %1=at least X% soft bottom, and above min depth

%filter_habitat_rows=study_area_polygons_PacCoastFisheryGIS_NUM(:,col_soft_hab_area)>0&study_area_polygons_PacCoastFisheryGIS_NUM(:,col_mean_depth)>halibut_max_depth;
study_area_polygons_PacCoastFisheryGIS_NUM_hab_depth=study_area_polygons_PacCoastFisheryGIS_NUM(filter_habitat_rows,:);
study_area_polygons_PacCoastFisheryGIS_TXT_noheader_hab_depth=study_area_polygons_PacCoastFisheryGIS_TXT_noheader(filter_habitat_rows,:);

%must first set blanks (NaNs) to zero for each habitat type
soft_tmp=study_area_polygons_PacCoastFisheryGIS_NUM_hab_depth(:,col_soft_hab_area); soft_tmp(isnan(soft_tmp)==1)=0;
mixed_tmp=study_area_polygons_PacCoastFisheryGIS_NUM_hab_depth(:,col_mixed_hab_area); mixed_tmp(isnan(mixed_tmp)==1)=0;
%then add the two habitat types to calculate total habitat
habitat_area_i_raw_km2=soft_tmp+mixed_tmp; % km2 habitat per patch within depth constraint
habitat_area_i_raw=habitat_area_i_raw_km2.*1e6; %m2 soft bottom per patch within depth constraint
depth_i=study_area_polygons_PacCoastFisheryGIS_NUM_hab_depth(:,col_mean_depth); %mean depth in each soft bottom patch within depth constraint
depth_i(depth_i>0)=0; %set above-sea level (ie positive depths) to zero (ie sea level)

numpatches=length(find(filter_habitat_rows==1)); %count of #rows=#patches in halibut model

lat_lon_habitat_patches=study_area_polygons_PacCoastFisheryGIS_NUM_hab_depth(:,[col_lat col_lon]); %matrix of lat and lon coordinates to each patch with soft habitat
target_fid_hab_depth=study_area_polygons_PacCoastFisheryGIS_NUM_hab_depth(:,col_target_fid);
blockid_for_each_soft_depth_patch=study_area_polygons_PacCoastFisheryGIS_NUM_hab_depth(:,block10_id); %CDFW block id associated with each habitat patch
distance_to_port_for_each_soft_depth_patch=study_area_polygons_PacCoastFisheryGIS_NUM_hab_depth(:,col_distance_to_port); %distance of each habitat patch to the nearest port
MPA_for_each_soft_depth_patch=strcmp('Y', study_area_polygons_PacCoastFisheryGIS_TXT_noheader_hab_depth(:,col_MPA));
PLATFORM_for_each_soft_depth_patch=strcmp('Y', study_area_polygons_PacCoastFisheryGIS_TXT_noheader_hab_depth(:,col_PLATFORM));
NAVIGATION_for_each_soft_depth_patch=strcmp('Y', study_area_polygons_PacCoastFisheryGIS_TXT_noheader_hab_depth(:,col_NAVIGATION));
ANCHORAGE_for_each_soft_depth_patch=strcmp('Y', study_area_polygons_PacCoastFisheryGIS_TXT_noheader_hab_depth(:,col_ANCHORAGE));
MILITARY_for_each_soft_depth_patch=strcmp('Y', study_area_polygons_PacCoastFisheryGIS_TXT_noheader_hab_depth(:,col_MILITARY));
OTHER_for_each_soft_depth_patch=strcmp('Y', study_area_polygons_PacCoastFisheryGIS_TXT_noheader_hab_depth(:,col_OTHER));
%Create vectord indicating if trawling and rec fishing is allowed
SQ_1trawling_0notrawling_for_each_soft_detph_patch=strcmp('Y',study_area_polygons_PacCoastFisheryGIS_TXT_noheader_hab_depth(:,col_Trawlable)); %set trawlable patches=1
SQ_1recfish_0recfish_for_each_soft_detph_patch=strcmp('yes',study_area_polygons_PacCoastFisheryGIS_TXT_noheader_hab_depth(:,col_rec_fish)); %set recreational fishable patches=1

%Areas where fishing currently not allowed
tmp=SQ_1trawling_0notrawling_for_each_soft_detph_patch+SQ_1recfish_0recfish_for_each_soft_detph_patch; % tmp=0 indicates no fishing in that patch, tmp=1 fishing by one group, tmp=2 fishing by both comm and rec
SQ_1fishable_0notfishable_for_each_soft_depth_patch=(tmp>=1); %0=fishing excluxed; 1=fishing allowed
SQ_1fishable_0notfishable_for_each_soft_depth_patch_ORIGINAL=SQ_1fishable_0notfishable_for_each_soft_depth_patch; %DON'T CHANGE

%vector of ones for each habitat type
hard1_other0_patch_vector=study_area_polygons_PacCoastFisheryGIS_NUM(:,col_hard_hab_area)>0;
mixed1_other0_patch_vector=study_area_polygons_PacCoastFisheryGIS_NUM(:,col_mixed_hab_area)>0;
soft1_other0_patch_vector=study_area_polygons_PacCoastFisheryGIS_NUM(:,col_soft_hab_area)>0;
unknown1_other0_patch_vector=study_area_polygons_PacCoastFisheryGIS_NUM(:,col_unknown_hab_area)>0;

%vector of ones for developable patches
% dev1_undev0_patch_vector=strcmp(study_area_polygons_PacCoastFisheryGIS_TXT_noheader(:,col_class), 'dev');
NoFish=strcmp(study_area_polygons_PacCoastFisheryGIS_TXT_noheader(:,col_class), 'NoFish');
All=strcmp(study_area_polygons_PacCoastFisheryGIS_TXT_noheader(:,col_class), 'All');
aqua_tmp=NoFish+All;
aqua_tmp(aqua_tmp>1)=1;
dev1_undev0_patch_vector=logical(aqua_tmp);

if runsetparams_1full==1
    %Coastline
    map_border_buffer_prop=0.001; %proportion shoreline map larger than edge patches
    disp(['Make shoreline map ',num2str(map_border_buffer_prop*100),'% larger than MSP domain'])
    %Lats are negative
    disp(['Lat map left = ',num2str(min(lat_lon_msp_domain(:,1))*(1+map_border_buffer_prop))])
    disp(['Lat map right = ',num2str(max(lat_lon_msp_domain(:,1))*(1-map_border_buffer_prop))])
    %Lons are positive
    disp(['Lon map upper = ',num2str(max(lat_lon_msp_domain(:,2))*(1+map_border_buffer_prop))])
    disp(['Lat map lower = ',num2str(min(lat_lon_msp_domain(:,2))*(1-map_border_buffer_prop))])
    disp('Plug the above values into <http://www.ngdc.noaa.gov/mgg_coastline/> to generate Matlab export of msp_domain.dat shoreline coordinates')
    %load shoreline coordinates
    msp_domain = textread('msp_domain.dat');
    %set mapping symbols and line dimensions (MAY NEED TO CHANGE THESE NUMBERS, depending on computer screen size)
    coastlinewidth=2; %thickness of coastline line on maps
    patchmarkersize=9; %size of patch squares on maps (set this so that the patches to "fit" together like a perfect grid)
    % patchmarkerlinewidth=.01; %width of the edges of the markers
    %Because the SCB coast does not cross the left border of the map, the
    %patches end up being the left-most point (not good for display). So,
    %create a point that is slightly left of the patches and have that be the edge :)
    lat_lon_leftmapedge=[min(lat_lon_msp_domain(:,1))*(1+map_border_buffer_prop) max(lat_lon_msp_domain(:,2))*(1+map_border_buffer_prop)];

        %Larval dispersal
    %ROMS-based connectivity matrix, given species (halibut) planktonic larval
    %dispersal (PLD) period, spawning season, and pattern of diel vertical migration
    %find the rows & cols corresponding to the soft bottom AND depth contraint (ie the rows indicated by target_fid_soft_depth)
    col_fid_labels=1; %column with patch number = fid value
     %new matrix with all cells in the domain (so, need to filter wrt habitat cells)
    load(strcat(input_data_dir,'halibut_connectivity_matrix_June27')) %6425x6425 cells
    cells_with_hab_NUM=study_area_polygons_PacCoastFisheryGIS_NUM(filter_habitat_rows,col_fid_labels); %fids of cells with suitable habitat

    %Empirical CDFW spatial block landings data
    load(strcat(input_data_dir,'PacCoastFisheryGIS_kglandings_by_block_Com_plus_Rec_sum2one'),'PacCoastFisheryGIS_kglandings_by_block_Com_plus_Rec_sum2one','unique_block_id','block_kglandings_inpatches_sum2one_Rec_Com')
    actual_comm_plus_rec_kglandings=177779; %Actual CDFW data on kg landings of CA halibut in SoCal (see Maunder et al. 2011 CA halibut stock assessment figures.xlsx sheet "Fig S2 + Fig S4 Other")
end

%demographics: Stock assessment refs refer to "Maunder et al. 2011 CA halibut Stock assessment <http://www.dfg.ca.gov/marine/sfmp/halibut-assessment.asp>"
%Growth
% Loo=137; %L infinity [cm], Rassweiller et al. 2012
% k=0.08; %growth rate [year^-1], Rassweiller et al. 2012
% t0=-1.2; % [year],relates to the size at first settlement, Rassweiller et al. 2012
%VBGF params (average of M and F) in Table A4.2.7. von Bertalanffy growth curve parameters in California Halibut Stock Assessment Background Information.pdf)
% Note: do NOT use data in Table A4.2.1. in California Halibut Stock Assessment Background Information.pdf, because not fitted to stock assessment model
Loo=126.25; %L infinity [cm],
k=0.1035; %growth rate [year^-1],
t0=-0.57; % [year],relates to the size at first settlement,
% C1=8.70E-06; % =C1, multiplier in length-weight equation [cm,kg], Rassweiller et al. 2012
% C2=3.05; % =C2, exponent in length-weight equation, Rassweiller et al. 2012
%L-W params (average of M and F) in Table B1.3.1. in Southern California Stock Assessment of California halibut.pdf and Table A4.2.1. in California Halibut Stock Assessment Background Information.pdf
C1=0.000008495; % =C1, multiplier in length-weight equation [cm,kg],
C2=3.03305; % =C2, exponent in length-weight equation,
%Life history
% max_age=30; %[years], Rassweiller et al. 2012r
max_age=27; %[years]; %max age (amax) is 30 and 23 for males and females (Table A4.9.1. in California Halibut Stock Assessment Background Information.pdf)
age_mature=4; %[years], Rassweiller et al. 2012
age_legal=5; %[years], Rassweiller et al. 2012; Legal harvest (commercial and recreational) is halibut at least 22 inches in length = 55 cm = 5 years given the stock assessment VBGF params
age_move=5; %Minimum age class at which individuals can move patches (not including dispersal during the larval stage).
% Domeier and Chun 1995: "California halibut larger than 500 mm total
% length (TL) tended to travel markedly greater distances
% than halibut smaller than 500 mm TL." Based on VBGF params from Maunder et al. 2011, 4 years = 47.58cmTL
% and 5 years=55.31cmTL. So, begin really moving when reach 5 years; age_move=5 years
%Density dependence
h_target=0.80; %median value for Family Pleuronidae (Myers 1999), refered in Maunder et al 2011
CRtarget= (4*h_target)/(1-h_target); %Hilborn and Walters 1992, White et al. 2010 Adapting the steepness parameter

%natural mortality
% delta=0.23.*ones(numpatches,1); %natural mortality rate [year^-1], Rassweiller et al. 2012
delta=0.25.*ones(numpatches,1); %natural mortality rate [year^-1]; %natural mortality (average of M and F). Natural mortality (M) is fixed at 0.2 for females and 0.3 for males based on assumptions made for summer flounder (Southern California Stock Assessment of California halibut.pdf)

% for every row (halibut habitat pathc) in target_fid_soft_depth
connect_matrix_rowcols_wrt_soft_depth=NaN(length(target_fid_hab_depth),1);
for index=1:length(target_fid_hab_depth)
    %find the row in Rachel's file that corresponds with the row in target_fid_soft_depth
     connect_matrix_rowcols_wrt_soft_depth(index)=find(cells_with_hab_NUM==target_fid_hab_depth(index)); %should be just one value
end
%Pull just those rows and save them as the larval dispersal kernel to use in the model
Dii=connect_matrix(connect_matrix_rowcols_wrt_soft_depth,connect_matrix_rowcols_wrt_soft_depth);
%Note: It is not a closed system (rows of Dii do not sum to one), thus
%some/many larvae are 'lost' out of the MSP study domain

%demographic functions
%For plotting use Lj and Wj fcns of age=0-max. But for the model use only fcn of age=1-max
%age-lenght-weight conversion
Lj_start_at_0age=Loo.*(1-exp(-k.*((0:max_age)-t0))); %cm length at each age class (von Bertalanffy function)
Wj_start_at_0age=C1.*(Lj_start_at_0age.^C2); %kg weight at each age class. Used for plotting
%For the model get rid of the age=0 value
Lj=Lj_start_at_0age;
Lj(1)=[];
Wj=Wj_start_at_0age;
Wj(1)=[];
Wij=repmat(Wj,numpatches,1); %weight at each age class, replicated across patch

%Economics
Theta_Density_prop=0.1; % TEST. default=0.1; %MR=MC break-even density for stock effect. 10% is an average value in Rassweiller et al. 2012
price_per_pound=4.8378; %Halibut price per pound landed. Calculated as ex-vessel
% value divided by landings of commercial halibut in southern California by year,
% averaged across years 2004-2011, the years when the halibut stock assessment
% indicates the catch and biomass to be stable.
% price=convmass(price_per_pound,'kg','lbm'); %price per kg halibut
lb_per_kg=2.20462;
price=price_per_pound*lb_per_kg;
