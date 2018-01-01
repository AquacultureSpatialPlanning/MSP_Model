% Filtering and seed excercises
clear
close all

%Scaled (0-1, except halibut, which is min-1) value of each sector wrt each policy
load EFPayoff_a_X_wrt_DM           
% EFPayoff_a_M_wrt_DM      1x279936            2239488  double  
% EFPayoff_a_F_wrt_DM      1x279936            2239488  double  
% EFPayoff_a_K_wrt_DM      1x279936            2239488  double  
% EFPayoff_a_H_wrt_DM      1x279936            2239488  double  
% EFPayoff_a_V_wrt_DM      1x279936            2239488  double 
% EFPayoff_a_B_wrt_DM      1x279936            2239488  double              
% EFPayoff_a_D_wrt_DM      1x279936            2239488  double  

load('Policy_i_a.mat')
% load TOA_data aMatrix
sectors={'Mussel','Finfish','Kelp','Halibut','Viewshed','Benthic','Disease'};
%% Filtering excercise:
close all
disp('Filtering excercise --------------------------------')
% aquaculture sector acheives >X_aqua% (5%) of its maximum possible value
X_aqua=0.05;
disp(['FILTER: each aqua sector acheives >',num2str(100*X_aqua),'% of its max value'])

    EFPayoff_a_M_wrt_DM_f01=(EFPayoff_a_M_wrt_DM<=X_aqua);%0=fits criteria
    disp(['EFPayoff_a_M_wrt_DM_f01 contains ',num2str(length(find(EFPayoff_a_M_wrt_DM_f01==0))),' policies'])

    EFPayoff_a_F_wrt_DM_f01=(EFPayoff_a_F_wrt_DM<=X_aqua);%0=fits criteria
    disp(['EFPayoff_a_F_wrt_DM_f01 contains ',num2str(length(find(EFPayoff_a_F_wrt_DM_f01==0))),' policies'])

    EFPayoff_a_K_wrt_DM_f01=(EFPayoff_a_K_wrt_DM<=X_aqua);%0=fits criteria
    disp(['EFPayoff_a_K_wrt_DM_f01 contains ',num2str(length(find(EFPayoff_a_K_wrt_DM_f01==0))),' policies'])

% no existing sectors are impacted by >X_existing% (5%)
X_existing=0.05;
X_existing_remainder=1-X_existing;
disp(['FILTER: no existing sectors are impacted by >',num2str(100*X_existing),'% of its max value'])

    EFPayoff_a_H_wrt_DM_f01=(EFPayoff_a_H_wrt_DM<X_existing_remainder);%0=fits criteria
    disp(['EFPayoff_a_H_wrt_DM_f01 contains ',num2str(length(find(EFPayoff_a_H_wrt_DM_f01==0))),' policies'])

    EFPayoff_a_V_wrt_DM_f01=(EFPayoff_a_V_wrt_DM<X_existing_remainder);%0=fits criteria
    % EFPayoff_a_V_wrt_DM_f01=(EFPayoff_a_V_wrt_DM<0.8);%0=fits criteria  <--NOTE CRITERIA
    disp(['EFPayoff_a_V_wrt_DM_f01 contains ',num2str(length(find(EFPayoff_a_V_wrt_DM_f01==0))),' policies'])

    EFPayoff_a_B_wrt_DM_f01=(EFPayoff_a_B_wrt_DM<X_existing_remainder);%0=fits criteria
    disp(['EFPayoff_a_B_wrt_DM_f01 contains ',num2str(length(find(EFPayoff_a_B_wrt_DM_f01==0))),' policies'])

    EFPayoff_a_D_wrt_DM_f01=(EFPayoff_a_D_wrt_DM<X_existing_remainder);%0=fits criteria
    disp(['EFPayoff_a_D_wrt_DM_f01 contains ',num2str(length(find(EFPayoff_a_D_wrt_DM_f01==0))),' policies'])

%filter wrt all criteria
EFPayoff_a_ALL_wrt_DM_f01=EFPayoff_a_M_wrt_DM_f01+EFPayoff_a_F_wrt_DM_f01+EFPayoff_a_K_wrt_DM_f01+...
    EFPayoff_a_H_wrt_DM_f01+EFPayoff_a_V_wrt_DM_f01+EFPayoff_a_B_wrt_DM_f01+EFPayoff_a_D_wrt_DM_f01;
disp(['EFPayoff_a_ALL_wrt_DM_f01 contains ',num2str(length(find(EFPayoff_a_ALL_wrt_DM_f01==0))),' policies'])

%Plot filter sector values
EFPayoff_a_X_wrt_DM_filter=[...
EFPayoff_a_M_wrt_DM(EFPayoff_a_ALL_wrt_DM_f01==0);...
EFPayoff_a_F_wrt_DM(EFPayoff_a_ALL_wrt_DM_f01==0);...
EFPayoff_a_K_wrt_DM(EFPayoff_a_ALL_wrt_DM_f01==0);...
EFPayoff_a_H_wrt_DM(EFPayoff_a_ALL_wrt_DM_f01==0);...
EFPayoff_a_V_wrt_DM(EFPayoff_a_ALL_wrt_DM_f01==0);...
EFPayoff_a_B_wrt_DM(EFPayoff_a_ALL_wrt_DM_f01==0);...
EFPayoff_a_D_wrt_DM(EFPayoff_a_ALL_wrt_DM_f01==0);...
];

%policies for each filtered plan:
Policy_i_a_ALL_wrt_DM_f01=Policy_i_a(:,EFPayoff_a_ALL_wrt_DM_f01==0);
save EFPayoff_a_ALL_wrt_DM_f01 EFPayoff_a_ALL_wrt_DM_f01
save Policy_i_a_ALL_wrt_DM_f01 Policy_i_a_ALL_wrt_DM_f01
save EFPayoff_a_X_wrt_DM_filter EFPayoff_a_X_wrt_DM_filter

figure
hold on
bar([1:7],EFPayoff_a_X_wrt_DM_filter)
xlabel('Sector')
set(gca,'XTick',1:length(sectors))
set(gca,'XTickLabel',sectors)
ylabel('Value')
title([num2str(length(find(EFPayoff_a_ALL_wrt_DM_f01==0))),' filtered policies: aqua >',num2str(100*X_aqua),'% max, existing impacted by <=',num2str(100*X_existing),'%'])

% ID the unique policies that meet the filter
Policy_i_a_filter_pi=find(EFPayoff_a_ALL_wrt_DM_f01==0);
Policy_i_a_filter=Policy_i_a(:,Policy_i_a_filter_pi);
Policy_i_a_filter_trans=Policy_i_a_filter';
[Policy_i_a_filter_unique,IA,IC] = unique(Policy_i_a_filter_trans,'rows','stable');
disp(['EFPayoff_a_ALL_wrt_DM_f01 contains ',num2str(length(Policy_i_a_filter_unique(:,1))),' UNIQUE policies'])

%% Seeds wrt bray curtis dissimilarity
braycurtis_d=NaN(length(Policy_i_a_filter_unique(:,1))); %shell
for i=1:length(Policy_i_a_filter_unique(:,1))
    for j=1:length(Policy_i_a_filter_unique(:,1))
        x=Policy_i_a_filter_unique(i,:)';
        y=Policy_i_a_filter_unique(j,:)';
        braycurtis_d(i,j)=braycurtis_crow_v1(x,y);
    end
end
figure
imagesc(braycurtis_d)
title('bray curtis dissimilarity index')

%Identify clusters
Z=linkage(braycurtis_d);
figure
[H T]=dendrogram(Z);
ylabel('Disimilarity')
xlabel('Plans')
%2-D non-metric multidimensional scaling
%the plan on the outer edge of each cluster is the obvious choice for user
%selection if the aim is to cover the maximum amount of variation (Linke et
%al. 2011)
Y=mdscale(braycurtis_d,2);
figure
for y=1:length(Y)
    plot(Y(y,1),Y(y,2),'w.')
    hold on
    text(Y(y,1),Y(y,2),num2str(y))
end
%Clusters from dendrogram
c1=[23 24 7 8 17 18];
c2=[9 13 14 19 25 20 26 21 15];
c3=[2 3 10];
c4=[1 12 16 22];
c5=[4 11 6 5];
%plot convex hull around mdscale points
for c=1:5
    eval(['ci=c',num2str(c),';'])
    x=Y(ci,1);
    y=Y(ci,2);  
    k = convhull(x,y);
    plot(x(k),y(k),'r-')
end
title('2-D NMDS of plans')
%Set axis the same for y and x
scaler=100;
% seed_plans_max_var=[2 6 22 23 26];
seed_plans_max_var=[6 2 26 23 22]; %rearrange order so bars in barplot ascned in order

% braycurtis_d_mean=mean(braycurtis_d);
% [B I]=sort(braycurtis_d_mean, 'descend')
% num_bc_set=5; % number of plans to plot
bc_set_u_i=seed_plans_max_var; %index among the 26 unique ones
bc_set_f_i=IA(bc_set_u_i); %index among the 26 unique of the 450 filtered ones
bc_set_policies=Policy_i_a_filter_trans(bc_set_f_i,:); %Actual policies
bc_set_policies=bc_set_policies'; %transpose so 1061x5
save bc_set_policies bc_set_policies
csvwrite('bc_set_policies.csv',bc_set_policies)
% bc_set_fullset_i=Policy_i_a_filter_pi(IC(bc_set_f_i)); %index among 279936 ones (not sure if this indexing is correct)
%Insert cell FIDs next to policies
load('Raw_Impacts_FID.mat') %structural array
Raw_Impacts
Raw_Impacts.FID;
csvwrite('FID1061cells.csv',Raw_Impacts.FID)
fid_bc_set_policies=[Raw_Impacts.FID bc_set_policies]; 
whos fid_bc_set_policies
save fid_bc_set_policies fid_bc_set_policies
csvwrite('fid_bc_set_policies.csv',fid_bc_set_policies)


% % close all
EFPayoff_a_X_wrt_DM_Seed_bc_set=[EFPayoff_a_X_wrt_DM_filter(:,bc_set_f_i)];%Actual payoffs
save EFPayoff_a_X_wrt_DM_Seed_bc_set EFPayoff_a_X_wrt_DM_Seed_bc_set

fig=figure
% subplot(2,1,2)
bar([1:7],[100.*EFPayoff_a_X_wrt_DM_Seed_bc_set])
% xlabel('Sector')
LEG=legend({'Seed 1','Seed 2','Seed 3','Seed 4','Seed 5'},'location','northwest','FontSize',24);
set(gca,'XTickLabel',sectors,'FontSize',18)
ylabel('Value [% of maximum]','FontSize',30)
% title('Seeds wrt bray-curtis')
set(gcf,'color','white'); 
box off
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
% fig.InvertHardcopy = 'off';
print(fig, 'Fig5a','-depsc','-tiff')

