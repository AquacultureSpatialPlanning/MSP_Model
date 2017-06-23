set(0,'DefaultTextFontSize',35)
set(0,'DefaultAxesFontSize',35)
set(0,'DefaultLineLineWidth',2)
set(0,'DefaultAxesLineWidth',2)
set(0,'DefaultSurfaceLineWidth',2)
grid on 
axis on

p = mfilename('fullpath');
part = fileparts(p);
parts = strsplit(part, '/Scripts');
DirPart = parts{end-1};
tic
disp('Evaluate MSP solutions')
addpath(genpath(DirPart))

load(strcat(DirPart,'/Input/Data/Raw_Impacts_FID'))
load(strcat(DirPart,'/Input/Data/tuned_params'))

view_MK = ((Raw_Impacts.Mussel > 0) + (Raw_Impacts.Kelp > 0)) > 0;
view_F = (Raw_Impacts.Finfish > 0);
MK_FID = Raw_Impacts.FID(view_MK);
F_FID = Raw_Impacts.FID(view_F);

tabulate(Raw_Impacts.Viewshed_Mussel_Kelp(view_MK) > 0)

MK_View_Fulldomain = Raw_Impacts.Viewshed_Mussel_Kelp(view_MK);
F_View_Fulldomain = Raw_Impacts.Viewshed_Finfish(view_F);

[c1,ia1,ib1] = intersect(target_fid_fulldomain,Raw_Impacts.FID(Raw_Impacts.Mussel > 0),'stable');
V1 = zeros(length(target_fid_fulldomain),1);V1(ia1) = 1;
[c1,ia1,ib1] = intersect(target_fid_fulldomain,Raw_Impacts.FID(Raw_Impacts.Finfish > 0),'stable');
V2 = zeros(length(target_fid_fulldomain),1);V2(ia1) = 1;
[c1,ia1,ib1] = intersect(target_fid_fulldomain,Raw_Impacts.FID(Raw_Impacts.Kelp > 0),'stable');
V3 = zeros(length(target_fid_fulldomain),1);V3(ia1) = 1;

%Viewshed_Mussel_Kelp = zeros(length(target_fid_fulldomain),1);Viewshed_Mussel_Kelp(V1 + V3 > 0) = Raw_Impacts.Viewshed_Mussel_Kelp;
%Viewshed_Finfish = zeros(length(target_fid_fulldomain),1);Viewshed_Finfish(V2 > 0) = Raw_Impacts.Viewshed_Finfish;
% Mussel/Kelp impact figure
    [c1,ia1,ib1] = intersect(target_fid_fulldomain,Raw_Impacts.FID(((Raw_Impacts.Mussel > 0) + (Raw_Impacts.Kelp > 0)) > 0));
    foo=zeros(length(target_fid_fulldomain),1);foo(ia1) = Raw_Impacts.Viewshed_Mussel_Kelp(ib1);
    patch_color=foo(ia1);
    tmp=lat_lon_msp_domain(ia1,:);
    figure
    scatter(tmp(:,1),tmp(:,2),patchmarkersize,patch_color,'s','filled') %all patches=grey filled in squares
    hold on
    scatter(lat_lon_leftmapedge(1),lat_lon_leftmapedge(2),0.1,'.','k')
    colorbar
    axis tight
%     xlabel('')
%     ylabel('')
%     title('Mussel and Kelp viewshed impacts')
    set(gcf,'color','white');
    plot(msp_domain(:,1),msp_domain(:,2),'k','linewidth',coastlinewidth) %coastline
    box on
    axis square
    axis tight
    print(strcat(DirPart,'/Input/Figures/S7A'),'-dpng')
% Finfish Impact Figures
    [c1,ia1,ib1] = intersect(target_fid_fulldomain,Raw_Impacts.FID(Raw_Impacts.Finfish > 0));
    foo=zeros(length(target_fid_fulldomain),1);foo(ia1) = Raw_Impacts.Viewshed_Finfish(ib1);
    patch_color=foo(ia1);
    tmp=lat_lon_msp_domain(ia1,:);
    figure
    scatter(tmp(:,1),tmp(:,2),patchmarkersize,patch_color,'s','filled') %all patches=grey filled in squares
    hold on
    scatter(lat_lon_leftmapedge(1),lat_lon_leftmapedge(2),0.1,'.','k')
    colorbar
    axis tight
%     xlabel('')
%     ylabel('')
%     title('')
    set(gcf,'color','white');
    plot(msp_domain(:,1),msp_domain(:,2),'k','linewidth',coastlinewidth) %coastline
    box on
    axis square
    axis tight
    print(strcat(DirPart,'/Input/Figures/S7B'),'-dpng')
