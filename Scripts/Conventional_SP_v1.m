%Calculates conventional management policies and outcomes

% from manuscript
% ratio = (aquaculture profit)/(scaled value of the most impacted existing sector) 
%from SI 7.2
% ratio to determine the suitability of each site for development 
% by each aquaculture sector: the annuity value of the aquaculture 
% sector if the site were developed, divided by the scaled value of 
% the most impacted existing sector at the site if that aquaculture 
% sector were developed there (Section 7.1). The rank order of sites 
% in relation to their suitability index was then used to simulate a 
% range of levels of development of aquaculture (1-1,061 of the developable 
% sites across the domain) under conventional planning. 

clear all
close all
tic %start timer

set(0,'defaultfigurecolor',[1 1 1])
set(0,'DefaultTextFontSize',10)
set(0,'DefaultAxesFontSize',10)
set(0,'DefaultLineLineWidth',2)
set(0,'DefaultAxesLineWidth',2)
set(0,'DefaultSurfaceLineWidth',2)

%% Calculate suitability index
%Value of each aqua in each site
%NUM1      1061x8             67904  double  
load('Raw_Impacts_FID')
% Raw_Impacts = 
%                      FID: [1061x1 int32]
%                   Mussel: [1061x1 double]
%                  Finfish: [1061x1 double]
%                     Kelp: [1061x1 double]
%                  Halibut: [1061x1 double]
%     Viewshed_Mussel_Kelp: [1061x1 double]
%         Viewshed_Finfish: [1061x1 double]
%          Benthic_Impacts: [1061x1 double]
%             Disease_Risk: [1061x1 double]
AquaValue=NaN(1061,8);
AquaValue(:,1)=Raw_Impacts.Mussel;
AquaValue(:,2)=Raw_Impacts.Finfish;
AquaValue(:,3)=Raw_Impacts.Kelp;

% scaled value of the most impacted existing sector at the site if that aquaculture sector were developed there
% X_n_i_p = Matrix of sector unitless value in each patch given policy choices
    % dimension 1 = i = 1....1061 sites
    % dimension 2 = n = 1....7 sectors
    % dimension 3 = p = 1...4 policies:
        %No development (p=1)
        %Mussel development (p=2)
        %Finfish development (p=3)
        %Kelp development (p=4)
% I = number of sites
load('TOA_data.mat', 'X_n_i_p','I')

SI=NaN(I,3); %Suitability index shell
%for each aqua sector
for a=1:3
    %for each site
    for i=1:I
        SI_numerator=AquaValue(i,a); %Aqua value at site
        X_existing_i_ND=X_n_i_p(i,3:7,1); %unitless value of existing sectors given No Development
        X_existing_i_aqua=X_n_i_p(i,3:7,a+1); %unitless value of existing sectors given aqua development a (note: a+1 because 1=ND)
        Ximpact_existing_i_aqua=X_existing_i_ND-X_existing_i_aqua; %Change in value (i.e, impact)
        SI_denominator=max(Ximpact_existing_i_aqua);
        SI(i,a)=SI_numerator/SI_denominator; %Note, if Aqua has not value at site, then SI=NaN
    end
end
%% Determine Conventional management plans
% Find max suitability values for each site and corresponding aqua type
[SI_Y,SI_I] = max(SI,[],2);
%Unconstrained
C_u=NaN(I,2); %shell of unconstrained conventional management policies. col 1= scenario #, col 2 = aqua policy
SI_Y_trimmer=SI_Y; %starts with all values; gets set to zero if site developed
for i=1:I %for each site developed (not necessarily in order of site numbers)
    i_aqua=find(SI_Y_trimmer==max(SI_Y_trimmer));
    i_aqua=i_aqua(1); %if multiple max (e.g., Inf) choose one
    SI_Y_trimmer(i_aqua)=0; %remove this site from the options to develop
    %Record the aqua type in the site location
    C_u(i_aqua,1)=i; %scenario #
    C_u(i_aqua,2)=SI_I(i_aqua)+1; %note +1 because SI is col=1-3 and aqua is M=2, F=3 and K=3
end
save C_u C_u
%Constrained - NOT CORRECT YET
C_c=NaN(I,2); %shell of unconstrained conventional management policies. col 1= scenario #, col 2 = aqua policy
SI_trimmer=SI;%starts with all values; row gets set to zero if site developed
index=0;
while index<I
% for i=1:I %for each site developed (not necessarily in order of site numbers)
%         SI_trimmer=SI;
        aqua_iter=ones(1,3); %keeps track of aqua chosen in each iteration (1=not chosen;0=chosen)
    for iter=1:3 %iterate through three aquas
        index=index+1; %
        if index==945
            945
%             pause(10)
        end
        aqua_iter_repmat=repmat(aqua_iter,I,1);
%         SI_search_iter=SI_search.*aqua_iter_repmat;
        SI_trimmer_iter=SI_trimmer.*aqua_iter_repmat;  
        [row col]=find(SI_trimmer_iter==max(SI_trimmer_iter(:)));
        row=row(1); col=col(1);
        SI_trimmer(row,:)=0; %remove this site from the options to develop
        aqua_iter(col)=0; %set aqua as chosen
        C_c(row,1)=index; %scenario #
        C_c(row,2)=col+1; %note +1 because SI is col=1-3 and aqua is M=2, F=3 and K=3
    end
end

%     
%        if last_policy==1; %ND
%         %choose any policy
%     elseif last_policy==2; %Mussel
%         %Set mussel suitability to zero
%         SI_search(:,1)=0;
%         SI_trimmer_search(:,1)=0;
%     elseif last_policy==3; %Finfish
%         %Set finfish suitability to zero
%         SI_search(:,2)=0;
%         SI_trimmer_search(:,2)=0;
%     elseif last_policy==4; %kelp
%         %Set kelp suitability to zero
%         SI_search(:,3)=0;
%         SI_trimmer_search(:,3)=0;
%     end        

