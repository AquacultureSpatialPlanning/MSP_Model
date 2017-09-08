%Calculate
close all
clear all
tic %start timer

% R=response (value or impact)
% V=value (positive numbers are beneficial to sector)
% X=scaled value (unitless, ranging 0-1, and sum to 1)
%
% n=sector (order: M, F, K, H, V, B, D)
% i=sites
% p=policy

%% load sector response data
% [NUM1,TXT1,RAW1]=xlsread('Raw_Patch_Data.csv'); %Joel says this data is rounded and thus INCORRECT
% 1 Mussel
% 2 Finfish
% 3 Kelp
% 4 Halibut
% 5 Views_Finfish
% 6 Views_Mussel_Kelp
% 7 Benthic
% 8 Disease

%INSTEAD, use .mat file provided by Joel:
load('~/MSP_Model/Scripts/CrowT0v1/Raw_Impacts.mat')
% Raw_Impacts =
%  1                 Mussel: [1061x1 double]
%  2                Finfish: [1061x1 double]
%  3                   Kelp: [1061x1 double]
%  4                Halibut: [1061x1 double]
%  5   Viewshed_Mussel_Kelp: [1061x1 double]
%  6       Viewshed_Finfish: [1061x1 double]
%  7        Benthic_Impacts: [1061x1 double]
%  8           Disease_Risk: [1061x1 double]
A = struct2cell(Raw_Impacts);
% Raw_Impacts converted into a matrix of doubles
%Note, NaN's already inserted for B and D (in 669 cells each)
%Note, order of V_M_K and V_F is reversed compared with Raw_Patch_Data.csv,
%so I changed the indexing in the below Response calculations
Raw_Impacts_double = cat(2,A{:});
NUM1=Raw_Impacts_double; %for use in below analysis

%% Parameters
N=7; %number of sectors
I=1061; %number of sites
P=4; %number of policies

%Calculate sector response in each patch to each policy
shell=zeros(I,N);
%No development (p=1)
R_n_i_p1=shell; %M, F, K, V, B and D zero (no value or no impact)
R_n_i_p1(:,4)=NUM1(:,4); %H full value
%Mussel development (p=2)
R_n_i_p2=shell; %F, K, H, B and D zero (no value or no impact)
R_n_i_p2(:,1)=NUM1(:,1); %M full value
% R_n_i_p2(:,5)=NUM1(:,6); %V partial impact (for analysis of Raw_Patch_Data.csv)
R_n_i_p2(:,5)=NUM1(:,5); %V partial impact (for analysis of Raw_Impacts.mat)
%Finfish development (p=3)
R_n_i_p3=shell; %M, K, and H zero (no value)
R_n_i_p3(:,2)=NUM1(:,2); %F full value
% R_n_i_p3(:,5)=NUM1(:,5); %V full impact (for analysis of Raw_Patch_Data.csv)
R_n_i_p3(:,5)=NUM1(:,6); %V full impact (for analysis of Raw_Impacts.mat)
R_n_i_p3(:,6)=NUM1(:,7); %B full impact
R_n_i_p3(:,7)=NUM1(:,8); %D full impact
%Kelp development (p=4)
R_n_i_p4=shell; %M, F, H, B and D zero (no value or no impact)
R_n_i_p4(:,3)=NUM1(:,3); %K full value
% R_n_i_p4(:,5)=NUM1(:,6); %V partial impact (for analysis of Raw_Patch_Data.csv)
R_n_i_p4(:,5)=NUM1(:,5); %V partial impact (for analysis of Raw_Impacts.mat)

%Create master matrix
R_n_i_p=NaN(I,N,P);
for p=1:P
    eval(['R_n_i_p(:,:,p)=R_n_i_p',num2str(p),';']);
end

%% Insert NaNs: this code necessary only for analysis of Raw_Patch_Data.csv
% % At sites where finfish cannot be developed, R=0 for B and D, indicating no impact
% % to B and D if F were developed there. This is incorrect. We didn't
% % calculate R at those sites, but it doesn't matter because they'll never
% % be developed. The issue can be corrected by simply setting R=NaN for B and D at
% % those sites, thereby preventing them from affecting the calculation of
% % V and X for B and D.
% ZeroF_i=find(NUM1(:,2)==0); %R_n_i_p3(:,2)==0 % Sites where F cannot be developed
% for p=1:P
%    R_n_i_p(ZeroF_i,6,p)=NaN; %Benthic
%    R_n_i_p(ZeroF_i,7,p)=NaN; %Disease
% end
%%

%Calculate maximum response by each
% sector across all sites and policies
tmp=nanmax(R_n_i_p,[],3); %first max across policies
R_bar_n=nanmax(tmp,[],1); %then max across sites

%Calculate values
V_n_i_p=NaN(I,N,P);
for n=1:N %for each sector
    for p=1:P %for each policy
        if n<=4 %3 M, F, K or H
            V_n_i_p(:,n,p)=R_n_i_p(:,n,p); %value
        elseif n>4 %V, B, D
            V_n_i_p(:,n,p)=R_bar_n(n) - R_n_i_p(:,n,p); %value
        end
    end
end

%Calculate unitless values
max_p_V_n_i_p=nanmax(V_n_i_p,[],3); %for each sector and site, the policy with the max value
sum_i_max_p_V_n_i_p=nansum(max_p_V_n_i_p,1); %sum of the above max values
X_n_i_p=NaN(I,N,P);
for n=1:N %for each sector
    for p=1:P %for each policy
        for i=1:I %for each site
            X_n_i_p(i,n,p)=V_n_i_p(i,n,p)./sum_i_max_p_V_n_i_p(n); %unitless value
        end
    end
end

%set up alpha matrix
a_range=linspace(0,1,6);
aMatrix=NaN(length(a_range)^N,N);
iMatrix=0;
aMi=0;
for aM=a_range
    aMi=aMi+1;
    disp(['aM=',num2str(aM),])
    for aF=a_range
        for aK=a_range
            for aH=a_range
                for aV=a_range
                    for aB=a_range
                        for aD=a_range
                            iMatrix=iMatrix+1;
                            aMatrix(iMatrix,:)=[aM,aF,aK,aH,aV,aB,aD];
                        end
                    end
                end
            end
        end
    end
end

%For each weighting scenario, determine optimal policy for each site
Policy_i_a=NaN(I,length(aMatrix));
iM_counter=round(linspace(1,length(aMatrix),10));
for iM=1:length(aMatrix)
    if isempty(intersect(iM,iM_counter))==0
        disp(['iM = ',num2str(iM)])
    end
    for i=1:I
       X_n_p_at_i=(squeeze(X_n_i_p(i,:,:)))'; %Unitless value of each policy to each sector
       aMatrix_repmatp=repmat(aMatrix(iM,:),P,1); %Weights of each sector, replicated across policies
       WX_n_p=X_n_p_at_i.*aMatrix_repmatp; %weighted unitless value of each policy to each sector
       sumWX_n_p=nansum(WX_n_p,2); %sum of weighted unitless values across sectors
       tmp=find(sumWX_n_p==max(sumWX_n_p)); %Identify policy that maximizes summed weighted unitless value
       Policy_i_a(i,iM)=tmp(1);
    end
end
save TOA_data.mat %save variables
save('Policy_i_a.mat','Policy_i_a','-v7.3');%save optimal policies (large file, so seperate)
disp(['Took ',num2str(toc/60/60),' hours']) %report run time

%% Analysis of results
figure
hist(aMatrix(:),6')
xlabel('Weight')
ylabel('Count')
title('All weighting scenarios')

figure
hist(Policy_i_a(:))
xlabel('Policy (1=no devel; 2=M; 3=F; 4=K)')
ylabel('Count')
title('All sites, all weighting scenarios')

figure
hist(R_n_i_p(:))
xlabel('Sector responses')
ylabel('Count')
title('All sites, all policies')

figure
hist(V_n_i_p(:))
xlabel('Sector values')
ylabel('Count')
title('All sites, all policies')

figure
hist(X_n_i_p(:))
xlabel('Sector unitless values')
ylabel('Count')
title('All sites, all policies')
