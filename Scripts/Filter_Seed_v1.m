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
load TOA_data aMatrix
%% Filtering excercise:
disp('Filtering excercise --------------------------------')
% aquaculture sector acheives >X_aqua% (5%)of its maximum possible value
X_aqua=0;
disp(['FILTER: each aqua sector acheives >',num2str(100*X_aqua),'% of its max value'])

    EFPayoff_a_M_wrt_DM_f=(EFPayoff_a_M_wrt_DM<=X_aqua);%0=fits criteria
    disp(['EFPayoff_a_M_wrt_DM_f contains ',num2str(length(find(EFPayoff_a_M_wrt_DM_f==0))),' policies'])

    EFPayoff_a_F_wrt_DM_f=(EFPayoff_a_F_wrt_DM<=X_aqua);%0=fits criteria
    disp(['EFPayoff_a_F_wrt_DM_f contains ',num2str(length(find(EFPayoff_a_F_wrt_DM_f==0))),' policies'])

    EFPayoff_a_K_wrt_DM_f=(EFPayoff_a_K_wrt_DM<=X_aqua);%0=fits criteria
    disp(['EFPayoff_a_K_wrt_DM_f contains ',num2str(length(find(EFPayoff_a_K_wrt_DM_f==0))),' policies'])

% no existing sectors are impacted by >X_existing% (5%)
X_existing=0.25;
X_existing_remainder=1-X_existing;
disp(['FILTER: no existing sectors are impacted by >',num2str(100*X_existing),'% of its max value'])

    EFPayoff_a_H_wrt_DM_f=(EFPayoff_a_H_wrt_DM<X_existing_remainder);%0=fits criteria
    disp(['EFPayoff_a_H_wrt_DM_f contains ',num2str(length(find(EFPayoff_a_H_wrt_DM_f==0))),' policies'])

    EFPayoff_a_V_wrt_DM_f=(EFPayoff_a_V_wrt_DM<X_existing_remainder);%0=fits criteria
    % EFPayoff_a_V_wrt_DM_f=(EFPayoff_a_V_wrt_DM<0.8);%0=fits criteria  <--NOTE CRITERIA
    disp(['EFPayoff_a_V_wrt_DM_f contains ',num2str(length(find(EFPayoff_a_V_wrt_DM_f==0))),' policies'])

    EFPayoff_a_B_wrt_DM_f=(EFPayoff_a_B_wrt_DM<X_existing_remainder);%0=fits criteria
    disp(['EFPayoff_a_B_wrt_DM_f contains ',num2str(length(find(EFPayoff_a_B_wrt_DM_f==0))),' policies'])

    EFPayoff_a_D_wrt_DM_f=(EFPayoff_a_D_wrt_DM<X_existing_remainder);%0=fits criteria
    disp(['EFPayoff_a_D_wrt_DM_f contains ',num2str(length(find(EFPayoff_a_D_wrt_DM_f==0))),' policies'])

%filter wrt all criteria
EFPayoff_a_ALL_wrt_DM_f=EFPayoff_a_M_wrt_DM_f+EFPayoff_a_F_wrt_DM_f+EFPayoff_a_K_wrt_DM_f+...
    EFPayoff_a_H_wrt_DM_f+EFPayoff_a_V_wrt_DM_f+EFPayoff_a_B_wrt_DM_f+EFPayoff_a_D_wrt_DM_f;
disp(['EFPayoff_a_ALL_wrt_DM_f contains ',num2str(length(find(EFPayoff_a_ALL_wrt_DM_f==0))),' policies'])

% 
% %Viewshed
% a=EFPayoff_a_V_wrt_DM(EFPayoff_a_ALL_wrt_DM_f==0);
% max(a) %at best V is impacted 1-0.8175=0.1825
% length(find(a>=0.95))% 12 plans result in a V impact no more than 20%

% ID the unique policies that meet the filter
Policy_i_a_filter=Policy_i_a(:,EFPayoff_a_ALL_wrt_DM_f==0);
Policy_i_a_filter_trans=Policy_i_a_filter';
[Policy_i_a_filter_unique,IA,IC] = unique(Policy_i_a_filter_trans,'rows','stable');
disp(['EFPayoff_a_ALL_wrt_DM_f contains ',num2str(length(Policy_i_a_filter_unique(:,1))),' UNIQUE policies'])
