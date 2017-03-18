close all
clear all
disp('------------------------')
%Nomenclature of subscripts to parameters used below
% I = patch = rows
% j = age = columns
% tic;
addpath(genpath('C:\Users\Joel\Desktop\Thesis YTD\Code'))% <- Add in a prompt for a directory
%% Settings
warning('off'); %no warning e.g. when divide by zero
%Set defaults for figures
set(0,'DefaultTextFontSize',35)
set(0,'DefaultAxesFontSize',35)
set(0,'DefaultLineLineWidth',6)
set(0,'DefaultAxesLineWidth',2)
set(0,'DefaultSurfaceLineWidth',2)
% Description of the model 
% Develop general functions (cost models etc)

%%% To do tasks: 1) Clean up tuner, 2) Get rid of unused variables in the
%%% entire model(s) 3) Make functions for generic tasks NPV etc. 

%% Load static sector data
prompt_sector_load='Calculate Sector Values(Y) or Load Values(N)? ';
str_set=input(prompt_sector_load,'s');
if strcmp(str_set,'Y')||strcmp(str_set,'y')
    % Load data 
        [NUM,TXT,RAW]=xlsread('SeaGrant_allFINALdata_complete_2015.xlsx');
        col_lat=2; %x-axis = latitude
        col_lon=3; %y-axis = longitude
        col_FID=1;
    % Time dimension of anaylsis 
        T=10; %NPV time horizon 
        discount_rate=0.05; %discount rate (d=0.05)
        discount_factor=1./((1+discount_rate).^[1:T]);
        discount_rate_iy_aqua=repmat(discount_factor,length(target_fid_fulldomain),1);
    % Aquaculture 
        % Mussel 
            disp('Parameterizing Mussel...') 
            mussel_upfrontcostcol=32; % Load all the mussel variables from the SeaGrant_allFINALdata_complete_2015.xlsx
            mussel_yieldkg_col=33;
            mussel_annualcost_col=34;
            mussel_annualprofit_col=35;
            mussel_10yrNPV_est_col=36;
            p_mussel=3.3; %Price in dollars per kg
            mussel_upfrontcost=NUM(:,mussel_upfrontcostcol); % Cost in the first year per patch
            mussel_yieldkg=NUM(:,mussel_yieldkg_col); % Yield per patch per year 
            mussel_annualcost=NUM(:,mussel_annualcost_col); % Annual cost per patch for all years following first year 
            mussel_rev=mussel_yieldkg.*p_mussel; % Yield of patch multipied by supplied price value 
            mussel_rev_mat=repmat(mussel_rev,1,T); % Make a revenue matrix by [# of patches T]
            mussel_year_1_total_cost=mussel_upfrontcost+mussel_annualcost; % The first year has to account for the start up costs plus the operating costs for that single year %MAY NEED TO EDIT THIS SINCE THEY MAY NOT PRODUCE ANYTHING IN THE FIRST YEAR 
            mussel_cost_mat=[mussel_year_1_total_cost repmat(mussel_annualcost,1,T-1)]; % Make a matrix of all of year 1 cost and a replication of all subsequent years within the range T 
            mussel_annual_profit_Piy=mussel_rev_mat-mussel_cost_mat; % Take the difference of the rev matrix and cost matrix to get annual profit per patch per year 
            mussel_NPViy=mussel_annual_profit_Piy.*discount_rate_iy_aqua; % Multiply the profit matrix by the discounted rate per year to get the NPV per patch per year 
            mussel_NPVi=sum(mussel_NPViy,2); % NPV per patch summed across all years 
            mussel_NPV=sum(mussel_NPVi,1); % Total NPV 
            % Scaling
            mussel_NPVi(mussel_NPVi<0)=0; %  Set negative values to zeros 
            mussel_NPVi_scaled_soln=mussel_NPVi./sum(mussel_NPVi); % Value scaled for measuring the RELATIVE VALUE OF SOLUTIONS
            mussel_NPVi_indices=find(mussel_NPVi>0); % Indices of patches that are postive value 
            V1=mussel_NPVi>0==1; % Converting zero and postive value to 0 and 1 binary 
            P1_full_scaled_obj=mussel_NPVi_scaled_soln; % Scaled for the objective function (Value relative to MIN and MAX impacts)  
        % Finfish 
            disp('Parameterizing FinFish...') 
            Fish_yieldkg_col=41;
            Fish_totalcost_col=42;
            Fish_annualrev_col=44;
            Fish_yieldkg=NUM(:,Fish_yieldkg_col); % Yield per patch per year 
            Fish_totalcost=NUM(:,Fish_totalcost_col);  % Annual cost per patch for all years 
            Fish_annualrev=NUM(:,Fish_annualrev_col); % Annual revenue per patch for all years 
            p_Fish=8; %Price of per kg, NOTE: The set price has a ton to do with what pacthes are set up to be profitable, changing this we estimate will have a dramtic effect on the number and spatial lay out of FISH developable patches, SENSITIVITY ANALYSIS 
            Fish_rev=Fish_yieldkg.*p_Fish; % Yield of patch multipied by supplied price value 
            Fish_rev_mat=repmat(Fish_rev,1,T);  % Make a revenue matrix by [# of patches T]
            Fish_cost_mat=[repmat(Fish_totalcost,1,T)]; % Make a replicated matrix with all values per patch of total yearly cost replicated across all years in the value T 
            Fish_annual_profit_Piy=Fish_rev_mat-Fish_cost_mat; % Take the difference of the rev matrix and cost matrix to get annual profit per patch per year 
            Fish_NPViy=Fish_annual_profit_Piy.*discount_rate_iy_aqua;   % Multiply the profit matrix by the discounted rate per year to get the NPV per patch per year 
            Fish_NPVi=sum(Fish_NPViy,2); % NPV per patch summed across all years 
            Fish_NPV=sum(Fish_NPVi,1); % Total NPV 
            % Scaling 
            Fish_NPVi(Fish_NPVi<0)=0; % Set negative values to zeros 
            Fish_NPVi_scaled_soln=Fish_NPVi./sum(Fish_NPVi);  % Value scaled for measuring the RELATIVE VALUE OF SOLUTIONS
            Fish_NPVi_indices=find(Fish_NPVi>0); % Indices of patches that are postive value 
            V2=Fish_NPVi>0==1; % Converting zero and postive value to 0 and 1 binary 
            P2_full_scaled_obj=Fish_NPVi_scaled_soln; % Scaled for the objective function (Value relative to MIN and MAX impacts)  
        % Kelp 
            disp('loading Kelp...')
            kelp_yieldkg_col=37;
            kelp_annualcost_col=38;
            kelp_upfrontcost_col=39;
            kelp_10yrNPV_est_col=40;
            p_kelp=3; %Price in dollars per kg 
            kelp_yieldkg=NUM(:,kelp_yieldkg_col);% Yield per patch per year 
            kelp_annualcost=NUM(:,kelp_annualcost_col); % Annual cost per patch for all years following first year 
            kelp_upfrontcost=NUM(:,kelp_upfrontcost_col);  % Cost in the first year per patch
            kelp_rev=kelp_yieldkg.*p_kelp; % Yield of patch multipied by supplied price value 
            kelp_rev_mat=repmat(kelp_rev,1,T);  % Make a revenue matrix by [# of patches T]
            kelp_year_1_total_cost=kelp_upfrontcost+kelp_annualcost; % The first year has to account for the start up costs plus the operating costs for that single year %MAY NEED TO EDIT THIS SINCE THEY MAY NOT PRODUCE ANYTHING IN THE FIRST YEAR
            kelp_cost_mat=[kelp_year_1_total_cost repmat(kelp_annualcost,1,T-1)]; % Make a matrix of all of year 1 cost and a replication of all subsequent years within the range T 
            kelp_annual_profit_Piy=kelp_rev_mat-kelp_cost_mat;  % Take the difference of the rev matrix and cost matrix to get annual profit per patch per year 
            kelp_NPViy=kelp_annual_profit_Piy.*discount_rate_iy_aqua;   % Multiply the profit matrix by the discounted rate per year to get the NPV per patch per year 
            kelp_NPVi=sum(kelp_NPViy,2); % NPV per patch summed across all years 
            kelp_NPV=sum(kelp_NPVi,1);  % Total NPV 
            % Sorting and scaling 
            kelp_NPVi(kelp_NPVi<0)=0;%  Set negative values to zeros 
            kelp_NPVi_scaled_soln=kelp_NPVi./sum(kelp_NPVi); % Value scaled for measuring the RELATIVE VALUE OF SOLUTIONS
            kelp_NPVi_indices=find(kelp_NPVi>0);  % Indices of patches that are postive value 
            V3=kelp_NPVi>0==1; % Converting zero and postive value to 0 and 1 binary 
            P3_full_scaled_obj=kelp_NPVi_scaled_soln; % Scaled for the objective function (Value relative to MIN and MAX impacts)  
        % Aquaculture parameters
            Aqua_dev_1Y_0N=V1+V2+V3>0==1; % Binary indicating any form of aquaculture 
            Aqua_dev_indices=find(Aqua_dev_1Y_0N==1);
        % Plot results 
            Aquaculture_Maps
    % Existing sectors 
        % Halibut 
            disp('Parameterizing Halibut...')
            % Have option to run the model dynamically 
%                 prompt_halibut='Run full halibut model(Y) or load variables(N)? ';
%                 str_halibut=input(prompt_halibut,'s');
%                 if strcmp(str_halibut,'Y')||strcmp(str_halibut,'y')
                    Halibut_tuner_free_params_v4
                    close all 
%                 else
%                     load('tuned_params')
%                 end 
            Halibut_10yrNPVi=Yi_fulldomain_NPV;
            Halibut_10yrNPVi(Halibut_10yrNPVi<0)=0; %Set negative values to zero 
            Halibut_max=sum(Halibut_10yrNPVi);
            Halibut_10yrNPV_scaled=Halibut_10yrNPVi./sum(Halibut_10yrNPVi); % Value scaled for measuring the relative value of the solutions
            Halibut_10yrNPV_indices=find(Halibut_10yrNPVi>0); % Indices of patches that are postive value
            P4_raw=Halibut_10yrNPVi; % RAW VALUE  
            V4=Halibut_10yrNPVi>0==1; % Converting zero and postive value to 0 and 1 binary 
            P4_full_devi=(Aqua_dev_1Y_0N==0).*P4_raw;% min raw Halibut value, due to full aqua development per patch 
            P4_full_dev=sum(P4_full_devi);% min raw Halibut value, due to full aqua development 
            P4_aqua_patches=Aqua_dev_1Y_0N.*P4_raw; %All halibut patches that can be impacted by aqua development 
            P4_aqua_patches_scaled_obj=P4_aqua_patches./sum(P4_aqua_patches); % Scaled for the objective function (Value relative to MIN and MAX impacts)  
            P4_max_dev_scaled_soln=P4_raw./sum(P4_raw); % Scaled for the solutions 
        % Viewshed 
            disp('Parameterizing Viewshed...');
            max_view_impact_aquadesign=zeros(length(target_fid_fulldomain),1);
            max_view_impact_aquadesign(V1==1)=1;
            max_view_impact_aquadesign(V3==1)=1; % May be overiding mussel cells, does not matter due to assignment of 1 to both 
            max_view_impact_aquadesign(V2==1)=2; % Overiding mussel/kelp with finfish where it can go 
            Viewshed_impacts_finfish_col_res=46;
            Viewshed_impacts_finfish_res=NUM(:,Viewshed_impacts_finfish_col_res); % Residential viewshed impacts for finfish 
            Viewshed_impacts_finfish_col_park=49;
            Viewshed_impacts_finfish_park=NUM(:,Viewshed_impacts_finfish_col_park); % Recreational viewshed impacts for finfish 
            Viewshed_impacts_kelpmussel_col_res=47;
            Viewshed_impacts_kelpmussel_res=NUM(:,Viewshed_impacts_kelpmussel_col_res); % Residential viewshed impacts for kelp/mussel
            Viewshed_impacts_kelpmussel_col_park=48;
            Viewshed_impacts_kelpmussel_park=NUM(:,Viewshed_impacts_kelpmussel_col_park);% Recreational viewshed impacts for kelp/mussel
            Viewshed_impacts_finfish=Viewshed_impacts_finfish_res+Viewshed_impacts_finfish_park;% Need to then sum the residential and recreational viewshed columns together for both kind of both kinds of impacted aquaculture 
            Viewshed_impacts_kelpmussel=Viewshed_impacts_kelpmussel_res+Viewshed_impacts_kelpmussel_park; %NOTE: Both respective viewshed vectors for mussel/kep and finfish are now the sum of both the residential values and recreational values 
            Viewshed_impacts_finfish(Viewshed_impacts_finfish<0)=0;
            Viewshed_impacts_kelpmussel(Viewshed_impacts_kelpmussel<0)=0;
            Viewshed_max_impact=max([Viewshed_impacts_finfish, Viewshed_impacts_kelpmussel],[],2); % Worst case scenario for viewshed impact 
            Viewshed_max_impact_scaled=Viewshed_max_impact/sum(Viewshed_max_impact);
            % Scaled the raw impacts wrt the worst case scenario impact 
            % Assuming that all patches with postive viewshed impact values can be developed for aquaculture
            Viewshed_impacts_finfish_scaled_obj=Viewshed_impacts_finfish./sum(Viewshed_max_impact); % Does not sum to one because mussel/kelp might be possible where finfish isn't 
            Viewshed_impacts_kelpmussel_scaled_obj=Viewshed_impacts_kelpmussel./sum(Viewshed_max_impact); % Does not sum to one because even if mussel and kelp development is maxed, impact is not as high as finfish 
            Viewshed_impacts_finfish_scaled_soln=Viewshed_impacts_finfish_scaled_obj; % Vector for scaling solutions
            Viewshed_impacts_kelpmussel_scaled_soln=Viewshed_impacts_kelpmussel_scaled_obj;% Vector for scaling solutions
            % Viewshed_impacts_indices=find(Viewshed_impacts>0);
            P5_wrt_finfish=Viewshed_impacts_finfish_scaled_obj; % Scaled for the objective function (Value relative to MIN and MAX impacts)  
            P5_wrt_kelpmussel=Viewshed_impacts_kelpmussel_scaled_obj; % Scaled for the objective function (Value relative to MIN and MAX impacts)  
            tmp=Viewshed_max_impact>0==1;
            V5=tmp;
        % Benthic Impacts 
            disp('Parameterizing Enviornmental Risk...');
            enviro_impact_col=50;
            enviro_impact=NUM(:,enviro_impact_col);
            fish_environ_logical=zeros(size(Fish_NPVi));
            fish_environ_logical(Fish_NPVi_indices)=1; % NOTE: Enviromental Impacts are only associated with patches that are developed for finfish aquaculture 
            enviro_impact(enviro_impact<0)=0;
            enviro_impact=enviro_impact.*fish_environ_logical;%NOTE: PREVIOUSLY PATCHES THAT HAD ENVIORN IMPACTS WERE INCLUDED EVEN IF THE VALUE OF FISH WAS NEGATIVE!!! 
            enviro_impact_scaled_soln=enviro_impact./sum(enviro_impact);
            % Assuming that all patches with environmental impacts postive values can be developed for finfish 
            P6_aqua_patches_scaled_obj=enviro_impact_scaled_soln;
            % P6=enviro_impact_scaled_soln;
        % Disease 
            max_impact_fish_X_logical=V2(Aqua_dev_indices);
            max_impact_fish_X=max_impact_fish_X_logical.*3;
            X=max_impact_fish_X;
            % Calculate Disease Max 
            X_domain=zeros(size(target_fid_fulldomain));
            X_domain(Aqua_dev_indices)=X';
            disease_vector_tmp=[(X_domain==3)+(X_domain==8)+(X_domain==13)+(X_domain==17)];
            disease_vector=disease_vector_tmp(Fish_NPVi_indices);
            tmp1=Fish_NPVi_indices(disease_vector==1); %of the whole domain, cells with aqua dev
            connect_matrix_disease=disease_connect_matrix(tmp1,tmp1);
            adj=connect_matrix_disease;
            eigen_disease=eigencentrality(adj); % Note that eignevalues are sometimes negative, so we need to convert it to postive 
            disease_metrici=abs(eigen_disease);
            disease_metric=sum(disease_metrici);
            disease_maxi=disease_metrici;
            disease_max=disease_metric; 
            % Calculate Disease Min 
            eigen_disease=eigencentrality(1); % Feeding a single patch for the minimum, note that eigen centraility will always be 1 for a single node network 
            disease_metrici=abs(eigen_disease);
            disease_metric=sum(disease_metrici);
            disease_mini=disease_metrici;
            disease_min=disease_metric; 
        % Plot results <- Need to edit the Halibut sections 
            Existing_Sector_Maps
else 
    load('Tuner_save_Jun_22_2015_15_28','mussel_NPVi_indices','Viewshed_impacts_kelpmussel_scaled_obj','Viewshed_impacts_finfish_scaled_obj','Viewshed_impacts_finfish_scaled_soln','Viewshed_impacts_kelpmussel_scaled_soln','disease_max','Yi_fulldomain_NPV','Halibut_max','P4_raw','Viewshed_impacts_finfish','Viewshed_impacts_kelpmussel','disease_connect_matrix','Fish_NPVi_indices','target_fid_fulldomain','Aqua_dev_indices','mussel_NPVi','Fish_NPVi','kelp_NPVi','Halibut_10yrNPVi','Viewshed_max_impact','enviro_impact','V2','disease_metrici')
end
%% Scale each sector 
%% Develop socioeconomic weighting preferences 
%% MSP Model 
    %this code show in concept how to conduct a 7D tradeoff analysis.
    %the variables are modeled after the SCB aquaculture tradeoff study, but
    % Set the epsilon level 
        prompt='Epsilon = ';
        str=input(prompt,'s');
        epsilon=str2double(str);
        N=1061; % Number of patches that can be developed for aqua 
        nS=7; % Number of sectors
    %Raw values of each sector in each patch
    % Emerging sectors
        M=mussel_NPVi(Aqua_dev_indices); % 90% of patches can be developed for M
        F=Fish_NPVi(Aqua_dev_indices); % 30% of the patches can be developed for F, some are outside M
        K=kelp_NPVi(Aqua_dev_indices); % 10% of the patches can be developed for K. F and K do not overlap
    % Existing sectors
        [~,IF,~] = intersect(Aqua_dev_indices,Fish_NPVi_indices); % Find out what indices of the 1061 developable cells are for finfish only 
        [~,IKM,~] = intersect(Aqua_dev_indices,mussel_NPVi_indices); % Do the same for kelp/mussel
        [~,IFKM,~] = intersect(Aqua_dev_indices(IF),Aqua_dev_indices(IKM));
        H = Halibut_10yrNPVi(Aqua_dev_indices); %Value that would be lost if M, F or K were developed in the patch
        V = Viewshed_max_impact(Aqua_dev_indices); % Max Impact 
        V_F=Viewshed_impacts_finfish; V_KM=Viewshed_impacts_kelpmussel; V_KM_diff=Viewshed_impacts_kelpmussel(Aqua_dev_indices);% Finfish and Kelp/Mussel Viewshed Impacts
        B = enviro_impact(Aqua_dev_indices); %Impact that would happen if F was developed there
        max_disease_adj=disease_connect_matrix(Fish_NPVi_indices,Fish_NPVi_indices);adj=max_disease_adj;eigen_max_disease=eigencentrality(adj);disease_metrici=abs(eigen_max_disease);
        D=zeros(length(Aqua_dev_indices),1);D(IF)=disease_metrici; %Impact that would happen if F was developed there
    %Matrix of the potential raw value of each sector in each patch
        EVr=[M, F, K, H, V, B, D];
        csvwrite('Raw_Patch_Data.csv',EVr);
    %Raw value of each sector in each patch wrt policy
        EVr1=EVr; EVr1(:,1:3)=0; %No development
        EVr2=EVr; EVr2(:,2:5)=0; %M develop
        EVr3=EVr; EVr3(:,1)=0; EVr3(:,3:7)=0; %F develop
        EVr4=EVr; EVr4(:,1:2)=0; EVr4(:,4:5)=0; %K develop
    %scale sectors so that their values remain proportional but sum to 1
        Ms=M./sum(M);
        Fs=F./sum(F);
        Ks=K./sum(K);
        Hs=H./sum(H);
        Vs=V./sum(V); Vsd=(V-V_KM_diff)/sum(V);
        Bs=B./sum(B); %first remove values that can't be affected
        Ds=D./sum(D);%first remove values that can't be affected
        Impacts=[Ms Fs Ks Hs Vs Bs Ds];
        csvwrite('All_Sectors_Impacts.csv',Impacts)
        clear mussel_NPVi_indices Viewshed_impacts_kelpmussel_scaled_obj Viewshed_impacts_finfish_scaled_obj Viewshed_impacts_finfish_scaled_soln Viewshed_impacts_kelpmussel_scaled_soln disease_max Yi_fulldomain_NPV Halibut_max P4_raw Viewshed_impacts_finfish Viewshed_impacts_kelpmussel disease_connect_matrix Fish_NPVi_indices target_fid_fulldomain mussel_NPVi Fish_NPVi kelp_NPVi Viewshed_max_impact enviro_impact V2 disease_metrici
    %Set sector weighting scenarios
    disp('Set sector weighting scenarios')
if epsilon==.2
        load('aMatrix')
        load('Dynamic_values')
else
        a_range=[0:epsilon:1]; %for test run use =[0 1]; for full run use=[0:0.2:1]. Or 0.1!!?? That would be amazing
        iMatrix=0;
        aMi=0;
        for aM=a_range
        aMi=aMi+1;
        disp(['aM=',num2str(aM),])
        aFi=0;
        for aF=a_range
        aFi=aFi+1;
        aKi=0;
        for aK=a_range
            aKi=aKi+1;
            aHi=0;
            for aH=a_range
                aHi=aHi+1;
                aVi=0;
                for aV=a_range
                    aVi=aVi+1;
                    aBi=0;
                    for aB=a_range
                        aBi=aBi+1;
                        aDi=0;
                        for aD=a_range
                            aDi=aDi+1;
                            iMatrix=iMatrix+1;
                            aMatrix(iMatrix,:)=[aM,aF,aK,aH,aV,aB,aD];
                        end
                    end
                end
            end
        end
        end
        end
        save('aMatrix','aMatrix')
    n_number_unique=unique_weights(aMatrix);
    disp(['Number of Unique Weights ',num2str(n_number_unique)])
    %Solve the sector weighting scenarios
    disp('Solve the sector weighting scenarios')
    shell=NaN(N,length(aMatrix(:,1)));
    Policy=shell;
%     EVwschoice=shell;
    % EVr_wrt_N_alpha=NaN(N,nS,length(aMatrix(:,1)));
    tic
    counter_display_numbers=round(linspace(1,length(aMatrix(:,1)),30));
    for ai=1:length(aMatrix(:,1)); 
    %Provide occasional update on progress
    %Note: this if..end command does slow the code speed some, but not much
    if ismember(ai,counter_display_numbers)==1
    disp(['ai=',num2str(ai),' of ',num2str(length(aMatrix(:,1))),' total scenarios done'])
    end
    %set particular scenario to be solved
    aM=aMatrix(ai,1);
    aF=aMatrix(ai,2);
    aK=aMatrix(ai,3);
    aV=aMatrix(ai,4);
    aH=aMatrix(ai,5);
    aB=aMatrix(ai,6);
    aD=aMatrix(ai,7);    

    %Calculate patch specific scaled and weighted values of the alternative policies
    EVsw1 = Hs.*aH + Vs.*aV + Bs.*aB + Ds.*aD; %No development
    EVsw2 = Ms.*aM + Vsd.*aV + Bs.*aB + Ds.*aD; %Develop M
    EVsw3 = Fs.*aF; %Develop F
    EVsw4 = Ks.*aK + Vsd.*aV + Bs.*aB + Ds.*aD; %Develop K

    %Identify which policy is best for each patch
    EVswoptions=[EVsw1 EVsw2 EVsw3 EVsw4];
    [Y, I] = max(EVswoptions,[],2); % I indicates which policy to implement in each patch, given the weights
    Policy(:,ai)=I-1; %rows=patch policy number; col=weighting scenario
    % EVwschoice(:,ai)=Y; %rows=patch policy weighted scaled value; col=weighting scenario
    end
    disp(['FINISHED: Took ',num2str(toc/60/60),' hours'])
    % Convert policy integers to ones which make more intuitive sense (i.e.
    % convert (1 to 0 for no development, etc)
    Static_plans=Policy;
%     Static_plans(Policy==1)=0;Static_plans(Policy==2)=1;Static_plans(Policy==3)=2;Static_plans(Policy==4)=3;
    Static_values=NaN(size(Static_plans,2),7);
    clear Policy
    % Find the values of each sector for each developed policy
    for index=1:size(Static_plans,2);
        X=zeros(6425,1);X(Aqua_dev_indices)=Static_plans(:,index);
        % Mussel 
            Static_values(index,1)=sum(Ms(Static_plans(:,index)==1));
        % Finfish
            Static_values(index,2)=sum(Fs(Static_plans(:,index)==2));
        % Kelp
            Static_values(index,3)=sum(Ks(Static_plans(:,index)==3));
        % Halibut 
            Static_values(index,4)=sum(Halibut_10yrNPVi(X==0))/sum(Halibut_10yrNPVi);
        % Viewshed
            Static_values(index,5)=1-(sum((V_F.*(X==2))+(V_KM.*[(X==1)+(X==3)]))/sum(V));
        % Benthic 
            Static_values(index,6)=sum(Bs(Static_plans(:,index)~=2));
        % Disease 
            Static_values(index,7)=sum(Ds(Static_plans(:,index)~=2));
    end
        save('Static_MSP_solutions','aMatrix','Static_plans','Static_values','-v7.3')
    %% Export Data to R 
        filename = 'C:\Users\Joel\Desktop\Thesis YTD\Code\MSP Planning Results April 2016';
        cd(filename)
        csvwrite('aMatrix.csv',aMatrix)
        csvwrite('Static_plans.csv',Static_plans)
        csvwrite('Static_values.csv',Static_values)
        aquaculture_maximum_alpha=0.05;
        existing_maxiumum_alpha=0.95;
        Case_study_data=casestudy_fx(aquaculture_maximum_alpha,existing_maxiumum_alpha,Static_plans,Static_values); % Sort for plans given inputted constraints on lines 164 and 165
        csvwrite('Static_plans_case_study.csv',Case_study_data.Static_plans_case_study)
        csvwrite('Static_values_case_study.csv',Case_study_data.Static_values_case_study)
        csvwrite('Static_percentage_case_study.csv',Case_study_data.Static_values_case_study)
        filename = 'C:\Users\Joel\Desktop\Thesis YTD\Code';
        cd(filename)
    %% Quality Control 
        figure;parallelcoords(Static_values,'Quantile',.5)
        figure;hist(Static_values(:,1));title('Mussl');figure;hist(Static_values(:,5));title('Finfish');figure;hist(Static_values(:,3));title('Kelp');
    %% Find Dynamic Solutions To the MSP solutions
    prompt_comp='Run Dynamic compiler(Y) or load date(N)? = ';
    str=input(prompt_comp,'s');
    if strcmp(str,'Y')
        % Number of Jobs 
            number_jobs=6;
        % Number of Tasks 
            number_tasks=6;
        % Number of 125
            number_iterations=7776;
        % Create Compiler 
            Static_compiler_matrix=NaN(number_jobs,number_tasks,number_iterations);
            itor_count1=0;
            itor_counter=0;
            while itor_count1<size(Static_compiler_matrix,1)
                itor_count1=itor_count1+1;
                itor_count2=0;
                while itor_count2<size(Static_compiler_matrix,2)
                    itor_count2=itor_count2+1;
                    itor_count3=0;
                    while itor_count3<size(Static_compiler_matrix,3)
                        itor_counter=itor_counter+1;
                        itor_count3=itor_count3+1;
                        Static_compiler_matrix(itor_count1,itor_count2,itor_count3)=itor_counter;
                    end
                end
            end 
    %         load('Static_MSP_solutions')
            disp('Starting Jobs......')
            c=parcluster('local');
            for index_jobs=1:size(Static_compiler_matrix,1);
                jfinish=createJob(c); 
                for index_tasks=1:size(Static_compiler_matrix,2)
                    index_iterations=1:size(Static_compiler_matrix,3);
                    timerVal = tic;
                    indexEF_tasks_compile=Static_compiler_matrix(index_jobs,index_tasks,index_iterations);
                    indexEF_tasks=indexEF_tasks_compile(:); X_in=Static_plans(:,indexEF_tasks);
                    field1='Task_solutions';value1=X_in;field2='Alphas';value2=aMatrix;field3='Job_index';value3=index_jobs;
                    field4='Task_index';value4=index_tasks;field5='Iteration_index';value5=index_iterations;field6='Index_mat';value6=Static_compiler_matrix;
                    tasks_in=struct(field1,value1,field2,value2,field3,value3,field4,value4,field5,value5,field6,value6);
                    t=createTask(jfinish,@Finishing_function_static, 1, {tasks_in}, 'CaptureDiary',true); 
                end
                jfinish.submit;
    %             disp(num2str(finishTime))
                disp(['Iteration set ',num2str(index_jobs),' submitted'])
            end
            disp('Job Finished...........')
             % Dynamic Job Compiler -> Add description later and finish later 
        % Must manually input the list 
        disp('Once jobs are finished uncomment and run the following section')
        disp('Right click each job in the scheduler and retrieve the outputs for each')
        disp('Update each of the job names in line 425')
        prompt_jobs='Enter Job ID use format [#s] = ';
        str_jobs=input(prompt_jobs,'s');
        epsilon=str2double(str);
        job_name_list={'job22_output','job23_output','job24_output','job25_output'};
        for index_job=1:number_jobs
            ID_tmp=ID_list(index_job);
            myvarname=sprintf('job%u_output',ID_tmp);
            myvarname=fetchOutputs(evalin('base',foo{1,1}));
            fetchOutputs(evalin('base',foo{1,1}))
        end
        Dynamic_values=NaN(number_jobs*number_tasks*number_iterations,29);
        counter=0;
        ID_list=[22,23,24,25];
        for index_job=1:number_jobs
            ID_tmp=ID_list(index_job);
            name_tmp={sprintf('job%u',ID_tmp)};
            outputs_tmp=fetchOutputs(evalin('base',name_tmp{1,1}));
            for index_task=1:number_tasks
                task_tmp=outputs_tmp{index_task,1}.Store;
                disp(['job# ',num2str(index_job),' task# ',num2str(index_task)])
                for index_iteri=1:number_iterations 
                    counter=counter+1;
                    for index_iterj=1:size(task_tmp,2)
                        Dynamic_values(counter,index_iterj)=task_tmp(index_iteri,index_iterj);
                    end
                end
            end
        end
        save('Dynamic_values','Dynamic_values')
    end
end
%% Conventionalv2 April 6th 2016
% Update to the conventional models, rewrote the sequential (now called
% constrained) model. Previous versions had empty cells in the higher
% iteration #'s 
    load('Tuner_save_Jun_22_2015_15_28','V1','V2','V3','target_fid_fulldomain','mussel_NPVi_indices','Viewshed_impacts_kelpmussel_scaled_obj','Viewshed_impacts_finfish_scaled_obj','Viewshed_impacts_finfish_scaled_soln','Viewshed_impacts_kelpmussel_scaled_soln','disease_max','Yi_fulldomain_NPV','Halibut_max','P4_raw','Viewshed_impacts_finfish','Viewshed_impacts_kelpmussel','disease_connect_matrix','Fish_NPVi_indices','target_fid_fulldomain','Aqua_dev_indices','mussel_NPVi','Fish_NPVi','kelp_NPVi','Halibut_10yrNPVi','Viewshed_max_impact','enviro_impact','V2','disease_metrici')
% Emerging sectors
    M=mussel_NPVi(Aqua_dev_indices); % 90% of patches can be developed for M
    F=Fish_NPVi(Aqua_dev_indices); % 30% of the patches can be developed for F, some are outside M
    K=kelp_NPVi(Aqua_dev_indices); % 10% of the patches can be developed for K. F and K do not overlap
% Existing sectors
    [~,IF,~] = intersect(Aqua_dev_indices,Fish_NPVi_indices); % Find out what indices of the 1061 developable cells are for finfish only 
    H = Halibut_10yrNPVi(Aqua_dev_indices); %Value that would be lost if M, F or K were developed in the patch
    V = Viewshed_max_impact(Aqua_dev_indices); % Max Impact 
    V_F=Viewshed_impacts_finfish; V_KM=Viewshed_impacts_kelpmussel;
    Vf_s=V_F(Aqua_dev_indices)./sum(V_F(Aqua_dev_indices));Vkm_s=V_KM(Aqua_dev_indices)./sum(V_KM(Aqua_dev_indices));
    B = enviro_impact(Aqua_dev_indices); %Impact that would happen if F was developed there
    max_disease_adj=disease_connect_matrix(Fish_NPVi_indices,Fish_NPVi_indices);adj=max_disease_adj;eigen_max_disease=eigencentrality(adj);disease_metrici=abs(eigen_max_disease);
    D=zeros(length(Aqua_dev_indices),1);D(IF)=disease_metrici; %Impact that would happen if F was developed there
%Matrix of the potential raw value of each sector in each patch
    Ms=M./sum(M);
    Fs=F./sum(F);
    Ks=K./sum(K);
    Hs=H./sum(H);
    Vs=V./sum(V); 
    Bs=B./sum(B); %first remove values that can't be affected
    Ds=D./sum(D);%first remove values that can't be affected
% Mussel impacted value 
    tmp=[Hs Vkm_s];[~,I_m] = max(tmp,[],2);Impact_M_tmp=NaN(length(tmp),1);
    for itor=1:length(tmp);Impact_M_tmp(itor)=tmp(itor,I_m(itor));end
    Impact_M=Impact_M_tmp;
    Mv_impacts=M./(Impact_M); Mv_impacts(isnan(Mv_impacts)==1)=0; [Y_M,I_M]=sort(Mv_impacts,'descend');
    Mv_impacts(isinf(Mv_impacts)==1)=M(isinf(Mv_impacts)==1);
% Finfish impacted value
    tmp=[Hs Vf_s Bs Ds];[~,I_f] = max(tmp,[],2);Impact_F_tmp=NaN(length(tmp),1);
    for itor=1:length(tmp);Impact_F_tmp(itor)=tmp(itor,I_f(itor));end
    Impact_F=Impact_F_tmp;
    Fv_impacts=F./(Impact_F); Fv_impacts(isnan(Fv_impacts)==1)=0; [Y_F,I_F]=sort(Fv_impacts,'descend');
    Fv_impacts(isinf(Fv_impacts)==1)=F(isinf(Fv_impacts)==1);
% Kelp impacted value
    tmp=[Hs Vkm_s];[~,I_k] = max(tmp,[],2);Impact_K_tmp=NaN(length(tmp),1);
    for itor=1:length(tmp);Impact_K_tmp(itor)=tmp(itor,I_k(itor));end
    Impact_K=Impact_K_tmp;
    Kv_impacts=K./(Impact_K); Kv_impacts(isnan(Kv_impacts)==1)=0; [Y_K,I_K]=sort(Kv_impacts,'descend');
    Kv_impacts(isinf(Kv_impacts)==1)=K(isinf(Kv_impacts)==1);
% Logical Vectors Representing the Subsets of the 1061 Developable Cells 
    % Which Can Be Developed for Each Aquaculture Farm Type 
    Vm=V1(Aqua_dev_indices);
    Vf=V2(Aqua_dev_indices);
    Vk=V3(Aqua_dev_indices);
% Use the cell FID as a common element between the farm types in the
% simulation 
    Aqua_FID=target_fid_fulldomain(V1+V2+V3>0);
    Aqua_FID_perm=Aqua_FID;
%% Unconstrained Model 
% Planning simulation loop 
    Suitability_Mat=[Mv_impacts Fv_impacts Kv_impacts];
    FID_Vector=Aqua_FID;
    FID_Vector_Reference=FID_Vector;
    Plan_initial=zeros(length(Aqua_dev_indices),1);
    Prev_plan=Plan_initial;
    tmp_plan=Prev_plan;
    itor=1;
    num_dev_patches=0;
    Unconstrained_Plan_Store=[Plan_initial NaN(length(Aqua_dev_indices))];
    Unconstrained_data=NaN(size(Unconstrained_Plan_Store,2),4);
    while num_dev_patches<length(Aqua_dev_indices)
        itor=itor+1;
        [Y,I]=sort(Suitability_Mat,'descend');
        [~,itor_aqua_type]=max(Y(1,:));
        Best_Cell_Y=Y(1,itor_aqua_type);
        Best_Cell_I=I(1,itor_aqua_type);
        Best_Cell_FID=FID_Vector(Best_Cell_I);
        [~,IA_1,~]=intersect(FID_Vector_Reference,Best_Cell_FID);% Identify the most suitable cell FID
        tmp_plan=Prev_plan;
        tmp_plan(IA_1)=itor_aqua_type;
        Unconstrained_Plan_Store(:,itor)=tmp_plan;
        [~,IA_2,~]=intersect(FID_Vector,Best_Cell_FID);
        FID_Vector(IA_2)=[];
        Suitability_Mat(IA_2,:)=[];
        num_dev_patches=sum(tmp_plan>0);
        disp(['iteration = ',num2str(itor)])
        disp(['Mussel = ',num2str(sum(tmp_plan==1)),' Finfish = ',num2str(sum(tmp_plan==2)),' Kelp = ',num2str(sum(tmp_plan==3))])
        if itor_aqua_type==1
         [~,IA_M_check,~]=intersect(target_fid_fulldomain(V1),FID_Vector_Reference(IA_1));
         if isempty(IA_M_check)==1
             disp('Error not developable cell for Mussels')
         end
        elseif itor_aqua_type==2
         [~,IA_F_check,~]=intersect(target_fid_fulldomain(V2),FID_Vector_Reference(IA_1));
         if isempty(IA_F_check)==1
             disp('Error not developable cell for Finfish')
         end  
        elseif itor_aqua_type==3
         [~,IA_K_check,~]=intersect(target_fid_fulldomain(V3),FID_Vector_Reference(IA_1));
         if isempty(IA_K_check)==1
             disp('Error not developable cell for Kelp')
         end
        end
        Unconstrained_data(itor,1)=itor;Unconstrained_data(itor,2)=sum(tmp_plan==1);Unconstrained_data(itor,3)=sum(tmp_plan==2);
        Unconstrained_data(itor,4)=sum(tmp_plan==3);
        Prev_plan=tmp_plan;
        disp(itor_aqua_type)
    end

        Unconstrained_Values_Static=NaN(size(Unconstrained_Plan_Store,2),7);
        for index=1:size(Unconstrained_Plan_Store,2);
            X=zeros(6425,1);X(Aqua_dev_indices)=Unconstrained_Plan_Store(:,index);
            % Mussel 
                Unconstrained_Values_Static(index,1)=sum(Ms(Unconstrained_Plan_Store(:,index)==1));
            % Finfish
                Unconstrained_Values_Static(index,2)=sum(Fs(Unconstrained_Plan_Store(:,index)==2));
            % Kelp
                Unconstrained_Values_Static(index,3)=sum(Ks(Unconstrained_Plan_Store(:,index)==3));
            % Halibut 
                Unconstrained_Values_Static(index,4)=sum(Halibut_10yrNPVi(X==0))/sum(Halibut_10yrNPVi);
            % Viewshed
                Unconstrained_Values_Static(index,5)=1-(sum((V_F.*(X==2))+(V_KM.*[(X==1)+(X==3)]))/sum(V));
            % Benthic 
                Unconstrained_Values_Static(index,6)=sum(Bs(Unconstrained_Plan_Store(:,index)~=2));
            % Disease 
                Unconstrained_Values_Static(index,7)=sum(Ds(Unconstrained_Plan_Store(:,index)~=2));
        end
%     prompt_comp='Run Dynamic compiler? = ';
%     str=input(prompt_comp,'s');
%     if strcmp(str,'Y')
        load('Tuner_save_condensed_April_8_2016')
        load('Tuner_save_Jun_22_2015_15_28','target_fid_fulldomain','Aqua_dev_indices','P1_full_scaled_obj','P2_full_scaled_obj','P3_full_scaled_obj','P4_aqua_patches_scaled_obj','Viewshed_impacts_finfish_scaled_obj',...
         'Viewshed_impacts_kelpmussel_scaled_obj','P6_aqua_patches_scaled_obj','disease_min','disease_max','Fish_NPVi_indices','SQ_1fishable_0notfishable_for_each_soft_depth_patch_ORIGINAL'...
         ,'Nij_initial_msy_TS','x0_PPUE_msy_TS','Fsum_msy','mussel_NPVi','Fish_NPVi','kelp_NPVi','P4_raw','Halibut_max','P4_full_dev','Viewshed_impacts_finfish','Viewshed_impacts_kelpmussel'...
         ,'enviro_impact','target_fid_hab_depth','Wij','numpatches','max_age','T','delta','age_legal','age_mature','Dii','alphaCR','beta_i','habitat_area_i','theta','price','gamma',...
        'distance_to_port_for_each_soft_depth_patch','Mii','habitat_area_i','age_move','Run_full_1_Run_dummy_0','dummy_indices','discount_rate_iy','disease_connect_matrix','x0_PPUE_initial','x0_PPUE_msy_TS')
        warning('off','all')
        gammatmp=gamma;
        Unconstrained_Dynamic_values=NaN(size(Unconstrained_Plan_Store,2),29);
        X_iter=Unconstrained_Plan_Store;
        for indexEF=1:length(Unconstrained_Dynamic_values)
        disp(['Unconstrained Iteration #',num2str(indexEF)])
%         X=zeros(6425,1);
%         X(Aqua_dev_indices)=X_iter(:,indexEF);
        X=X_iter(:,indexEF);
        X_domain=zeros(size(target_fid_fulldomain));
        X_domain(Aqua_dev_indices)=X;
        [P1,P2,P3,P4,P5,P6,P7,P1i,P2i,P3i,P4i,P5i,P6i,P7i,mussel_NPVi_tmp,mussel_NPV_tmp,Fish_NPVi_tmp,Fish_NPV_tmp,kelp_NPVi_tmp,kelp_NPV_tmp,Halibut_NPVi_tmp,Halibut_NPV_tmp,Halibut_NPV_scaled_tmp,Halibut_Yi_tmp,halibut_dynamic_scaled_tmp,Viewshed_raw_valuei_fin_tmp,Viewshed_raw_valuei_kelpmussel_tmp,Viewshed_raw_value_fin_tmp,Viewshed_raw_value_kelpmussel_tmp,enviro_impacti_tmp,enviro_impact_tmp,disease_rawi,disease_raw_sum]=Sector_calcs_wrt_X(X,target_fid_fulldomain,Aqua_dev_indices,P1_full_scaled_obj,P2_full_scaled_obj,P3_full_scaled_obj,P4_aqua_patches_scaled_obj,Viewshed_impacts_finfish_scaled_obj,...
        Viewshed_impacts_kelpmussel_scaled_obj,P6_aqua_patches_scaled_obj,disease_min,disease_max,Fish_NPVi_indices,SQ_1fishable_0notfishable_for_each_soft_depth_patch_ORIGINAL,Nij_initial_msy_TS,x0_PPUE_msy_TS,Fsum_msy,mussel_NPVi,Fish_NPVi,kelp_NPVi,P4_raw,Halibut_max,P4_full_dev,Viewshed_impacts_finfish,Viewshed_impacts_kelpmussel,enviro_impact,target_fid_hab_depth,Wij,numpatches,max_age,T,delta,age_legal,age_mature,Dii,alphaCR,beta_i,habitat_area_i,theta,price,gammatmp,distance_to_port_for_each_soft_depth_patch,Mii,age_move,Run_full_1_Run_dummy_0,dummy_indices,discount_rate_iy,disease_connect_matrix);
        Unconstrained_Dynamic_values(indexEF,1)=NaN;
        Unconstrained_Dynamic_values(indexEF,2)=NaN;
        Unconstrained_Dynamic_values(indexEF,3)=NaN;
        Unconstrained_Dynamic_values(indexEF,4)=NaN;
        Unconstrained_Dynamic_values(indexEF,5)=NaN;
        Unconstrained_Dynamic_values(indexEF,6)=NaN;
        Unconstrained_Dynamic_values(indexEF,7)=NaN;
        %% GA Outputs
        Unconstrained_Dynamic_values(indexEF,8)=NaN;
        Unconstrained_Dynamic_values(indexEF,9)=NaN;
        Unconstrained_Dynamic_values(indexEF,10)=NaN;
        Unconstrained_Dynamic_values(indexEF,11)=NaN;
        %% Scaled Sector Totals (Static)
        Unconstrained_Dynamic_values(indexEF,12)=P1; 
        Unconstrained_Dynamic_values(indexEF,13)=P2;
        Unconstrained_Dynamic_values(indexEF,14)=P3;
        Unconstrained_Dynamic_values(indexEF,15)=P4;
        Viewshed_rawi=(Viewshed_impacts_finfish.*(X_domain==2))+...
        (Viewshed_impacts_kelpmussel.*(X_domain==1)+(X_domain==3));
        Viewshed_Scaled_Impacts=sum(Viewshed_rawi)/sum(Viewshed_max_impact);
        Viewshed_Scaled=1-Viewshed_Scaled_Impacts;
        Unconstrained_Dynamic_values(indexEF,16)=Viewshed_Scaled; 
        Unconstrained_Dynamic_values(indexEF,17)=P6;
        Unconstrained_Dynamic_values(indexEF,18)=P7; 
        %% Sector Raw Totals 
        Unconstrained_Dynamic_values(indexEF,19)=mussel_NPV_tmp; % Mussel Raw 
        Unconstrained_Dynamic_values(indexEF,20)=Fish_NPV_tmp; % Finfish Raw 
        Unconstrained_Dynamic_values(indexEF,21)=kelp_NPV_tmp; % Kelp Raw 
        Unconstrained_Dynamic_values(indexEF,22)=Halibut_NPV_tmp; % Halibut Static Raw
        Unconstrained_Dynamic_values(indexEF,23)=Halibut_NPV_scaled_tmp; % Halibut Static Scaled 
        Unconstrained_Dynamic_values(indexEF,24)=Halibut_Yi_tmp; % Halibut Dynamic Raw 
        Unconstrained_Dynamic_values(indexEF,25)=Halibut_Yi_tmp/Halibut_max;% Halibut Dynamic Scaled 
        Unconstrained_Dynamic_values(indexEF,26)=Viewshed_raw_value_fin_tmp; % Viewshed Finfish Raw 
        Unconstrained_Dynamic_values(indexEF,27)=Viewshed_raw_value_kelpmussel_tmp; % Viewshed Kelp/Mussel Raw 
        Unconstrained_Dynamic_values(indexEF,28)=enviro_impact_tmp; % Enviormental Impacts Raw 
        Unconstrained_Dynamic_values(indexEF,29)=disease_raw_sum;% Disease Raw  
        end
        save('Unconstrained_Dynamic_Values','Unconstrained_Dynamic_values')
        csvwrite('Unconstrained_Dynamic_Values_April.csv',Unconstrained_Dynamic_values)
%     end
%% Constrained Model 
    Suitability_Mat=[Mv_impacts Fv_impacts Kv_impacts];
    FID_Vector=Aqua_FID;
    FID_Vector_Reference=FID_Vector;
    Plan_initial=zeros(length(Aqua_dev_indices),1);
    Prev_plan=Plan_initial;
    tmp_plan=Prev_plan;
    itor=1;
    num_dev_patches=0;
    Constrained_Plan_Store=NaN(length(Aqua_dev_indices));
    Constrained_data=NaN(size(Constrained_Plan_Store,2),4);
    while num_dev_patches<length(Aqua_dev_indices)
        n_aqua_types=find(sum(Suitability_Mat>0,1)>0);
        for index=1:length(n_aqua_types)
            itor_aqua_type=n_aqua_types(index);
            itor=itor+1;
            [Y,I]=sort(Suitability_Mat,'descend');
            Best_Cell_Y=Y(1,itor_aqua_type);
            Best_Cell_I=I(1,itor_aqua_type);
            Best_Cell_FID=FID_Vector(Best_Cell_I);
            [~,IA_1,~]=intersect(FID_Vector_Reference,Best_Cell_FID);% Identify the most suitable cell FID
            tmp_plan=Prev_plan;
            tmp_plan(IA_1)=itor_aqua_type;
            Constrained_Plan_Store(:,itor)=tmp_plan;
            [~,IA_2,~]=intersect(FID_Vector,Best_Cell_FID);
            FID_Vector(IA_2)=[];
            Suitability_Mat(IA_2,:)=[];
            num_dev_patches=sum(tmp_plan>0);
            disp(['iteration = ',num2str(itor)])
            disp(['Mussel = ',num2str(sum(tmp_plan==1)),' Finfish = ',num2str(sum(tmp_plan==2)),' Kelp = ',num2str(sum(tmp_plan==3))])
            if itor_aqua_type==1
             [~,IA_M_check,~]=intersect(target_fid_fulldomain(V1),FID_Vector_Reference(IA_1));
             if isempty(IA_M_check)==1
                 disp('Error not developable cell for Mussels')
             end
            elseif itor_aqua_type==2
             [~,IA_F_check,~]=intersect(target_fid_fulldomain(V2),FID_Vector_Reference(IA_1));
             if isempty(IA_F_check)==1
                 disp('Error not developable cell for Finfish')
             end  
            elseif itor_aqua_type==3
             [~,IA_K_check,~]=intersect(target_fid_fulldomain(V3),FID_Vector_Reference(IA_1));
             if isempty(IA_K_check)==1
                 disp('Error not developable cell for Kelp')
             end
            end
            Constrained_data(itor,1)=itor;Constrained_data(itor,2)=sum(tmp_plan==1);Constrained_data(itor,3)=sum(tmp_plan==2);
            Constrained_data(itor,4)=sum(tmp_plan==3);
            Prev_plan=tmp_plan;
            disp(itor_aqua_type)
        end
    end
    Constrained_Plan_Store(:,1)=zeros(length(Aqua_dev_indices),1);
    Constrained_Values_Static=NaN(size(Constrained_Plan_Store,2),7);
    for index=1:size(Constrained_Plan_Store,2);
        X=zeros(6425,1);X(Aqua_dev_indices)=Constrained_Plan_Store(:,index);
        % Mussel 
            Constrained_Values_Static(index,1)=sum(Ms(Constrained_Plan_Store(:,index)==1));
        % Finfish
            Constrained_Values_Static(index,2)=sum(Fs(Constrained_Plan_Store(:,index)==2));
        % Kelp
            Constrained_Values_Static(index,3)=sum(Ks(Constrained_Plan_Store(:,index)==3));
        % Halibut 
            Constrained_Values_Static(index,4)=sum(Halibut_10yrNPVi(X==0))/sum(Halibut_10yrNPVi);
        % Viewshed
            Constrained_Values_Static(index,5)=1-(sum((V_F.*(X==2))+(V_KM.*[(X==1)+(X==3)]))/sum(V));
        % Benthic 
            Constrained_Values_Static(index,6)=sum(Bs(Constrained_Plan_Store(:,index)~=2));
        % Disease 
            Constrained_Values_Static(index,7)=sum(Ds(Constrained_Plan_Store(:,index)~=2));
    end
%     prompt_comp='Run Dynamic compiler? = ';
%     str=input(prompt_comp,'s');
%     if strcmp(str,'Y')
        load('Tuner_save_condensed_April_8_2016')
        load('Tuner_save_Jun_22_2015_15_28','target_fid_fulldomain','Aqua_dev_indices','P1_full_scaled_obj','P2_full_scaled_obj','P3_full_scaled_obj','P4_aqua_patches_scaled_obj','Viewshed_impacts_finfish_scaled_obj',...
        'Viewshed_impacts_kelpmussel_scaled_obj','P6_aqua_patches_scaled_obj','disease_min','disease_max','Fish_NPVi_indices','SQ_1fishable_0notfishable_for_each_soft_depth_patch_ORIGINAL'...
        ,'Nij_initial_msy_TS','x0_PPUE_msy_TS','Fsum_msy','mussel_NPVi','Fish_NPVi','kelp_NPVi','P4_raw','Halibut_max','P4_full_dev','Viewshed_impacts_finfish','Viewshed_impacts_kelpmussel'...
        ,'enviro_impact','target_fid_hab_depth','Wij','numpatches','max_age','T','delta','age_legal','age_mature','Dii','alphaCR','beta_i','habitat_area_i','theta','price','gamma',...
        'distance_to_port_for_each_soft_depth_patch','Mii','habitat_area_i','age_move','Run_full_1_Run_dummy_0','dummy_indices','discount_rate_iy','disease_connect_matrix','x0_PPUE_initial','x0_PPUE_msy_TS')
        warning('off','all')
        gammatmp=gamma;
        Constrained_Dynamic_values=NaN(size(Constrained_Plan_Store,2),29);
        X_iter=Constrained_Plan_Store;
        for indexEF=1:length(Constrained_Dynamic_values)
        disp(['Constrained Iteration #',num2str(indexEF)])
        X=X_iter(:,indexEF);
        X_domain=zeros(size(target_fid_fulldomain));
        X_domain(Aqua_dev_indices)=X;
        [P1,P2,P3,P4,P5,P6,P7,P1i,P2i,P3i,P4i,P5i,P6i,P7i,mussel_NPVi_tmp,mussel_NPV_tmp,Fish_NPVi_tmp,Fish_NPV_tmp,kelp_NPVi_tmp,kelp_NPV_tmp,Halibut_NPVi_tmp,Halibut_NPV_tmp,Halibut_NPV_scaled_tmp,Halibut_Yi_tmp,halibut_dynamic_scaled_tmp,Viewshed_raw_valuei_fin_tmp,Viewshed_raw_valuei_kelpmussel_tmp,Viewshed_raw_value_fin_tmp,Viewshed_raw_value_kelpmussel_tmp,enviro_impacti_tmp,enviro_impact_tmp,disease_rawi,disease_raw_sum]=Sector_calcs_wrt_X(X,target_fid_fulldomain,Aqua_dev_indices,P1_full_scaled_obj,P2_full_scaled_obj,P3_full_scaled_obj,P4_aqua_patches_scaled_obj,Viewshed_impacts_finfish_scaled_obj,...
        Viewshed_impacts_kelpmussel_scaled_obj,P6_aqua_patches_scaled_obj,disease_min,disease_max,Fish_NPVi_indices,SQ_1fishable_0notfishable_for_each_soft_depth_patch_ORIGINAL,Nij_initial_msy_TS,x0_PPUE_msy_TS,Fsum_msy,mussel_NPVi,Fish_NPVi,kelp_NPVi,P4_raw,Halibut_max,P4_full_dev,Viewshed_impacts_finfish,Viewshed_impacts_kelpmussel,enviro_impact,target_fid_hab_depth,Wij,numpatches,max_age,T,delta,age_legal,age_mature,Dii,alphaCR,beta_i,habitat_area_i,theta,price,gammatmp,distance_to_port_for_each_soft_depth_patch,Mii,age_move,Run_full_1_Run_dummy_0,dummy_indices,discount_rate_iy,disease_connect_matrix);
        Constrained_Dynamic_values(indexEF,1)=NaN;
        Constrained_Dynamic_values(indexEF,2)=NaN;
        Constrained_Dynamic_values(indexEF,3)=NaN;
        Constrained_Dynamic_values(indexEF,4)=NaN;
        Constrained_Dynamic_values(indexEF,5)=NaN;
        Constrained_Dynamic_values(indexEF,6)=NaN;
        Constrained_Dynamic_values(indexEF,7)=NaN;
        %% GA Outputs
        Constrained_Dynamic_values(indexEF,8)=NaN;
        Constrained_Dynamic_values(indexEF,9)=NaN;
        Constrained_Dynamic_values(indexEF,10)=NaN;
        Constrained_Dynamic_values(indexEF,11)=NaN;
        %% Scaled Sector Totals (Static)
        Constrained_Dynamic_values(indexEF,12)=P1; 
        Constrained_Dynamic_values(indexEF,13)=P2;
        Constrained_Dynamic_values(indexEF,14)=P3;
        Constrained_Dynamic_values(indexEF,15)=P4;
        Viewshed_rawi=(Viewshed_impacts_finfish.*(X_domain==2))+...
        (Viewshed_impacts_kelpmussel.*(X_domain==1)+(X_domain==3));
        Viewshed_Scaled_Impacts=sum(Viewshed_rawi)/sum(Viewshed_max_impact);
        Viewshed_Scaled=1-Viewshed_Scaled_Impacts;
        Unconstrained_Dynamic_values(indexEF,16)=Viewshed_Scaled; 
        Constrained_Dynamic_values(indexEF,17)=P6;
        Constrained_Dynamic_values(indexEF,18)=P7; 
        %% Sector Raw Totals 
        Constrained_Dynamic_values(indexEF,19)=mussel_NPV_tmp; % Mussel Raw 
        Constrained_Dynamic_values(indexEF,20)=Fish_NPV_tmp; % Finfish Raw 
        Constrained_Dynamic_values(indexEF,21)=kelp_NPV_tmp; % Kelp Raw 
        Constrained_Dynamic_values(indexEF,22)=Halibut_NPV_tmp; % Halibut Static Raw
        Constrained_Dynamic_values(indexEF,23)=Halibut_NPV_scaled_tmp; % Halibut Static Scaled 
        Constrained_Dynamic_values(indexEF,24)=Halibut_Yi_tmp; % Halibut Dynamic Raw 
        Constrained_Dynamic_values(indexEF,25)=Halibut_Yi_tmp/Halibut_max; % Halibut Dynamic Scaled 
        Constrained_Dynamic_values(indexEF,26)=Viewshed_raw_value_fin_tmp; % Viewshed Finfish Raw 
        Constrained_Dynamic_values(indexEF,27)=Viewshed_raw_value_kelpmussel_tmp; % Viewshed Kelp/Mussel Raw 
        Constrained_Dynamic_values(indexEF,28)=enviro_impact_tmp; % Enviormental Impacts Raw 
        Constrained_Dynamic_values(indexEF,29)=disease_raw_sum;% Disease Raw  
        end
        save('Constrained_Dynamic_Values','Constrained_Dynamic_values')
        csvwrite('Constrained_Dynamic_Values_April.csv',Constrained_Dynamic_values)
%     end
    % Save variables
        filename = 'C:\Users\Joel\Desktop\Thesis YTD\Code\MSP Planning Results April 2016';
        cd(filename)
        save('Conventional_Results_April','Unconstrained_Plan_Store','Unconstrained_Values_Static','Constrained_Plan_Store','Constrained_Values_Static')
        csvwrite('Unconstrained_Plans_April.csv',Unconstrained_Plan_Store)
        csvwrite('Unconstrained_Static_Values_April.csv',Unconstrained_Values_Static)
        csvwrite('Constrained_Plans_April.csv',Constrained_Plan_Store)
        csvwrite('Constrained_Static_Values_April.csv',Constrained_Values_Static)
        filename = 'C:\Users\Joel\Desktop\Thesis YTD\Code';
        cd(filename)
    % Plot conventional management figures 
        % Unconstrained 
            figure
            plot(Unconstrained_data(:,1),Unconstrained_data(:,2),'-b')
            hold on 
            plot(Unconstrained_data(:,1),Unconstrained_data(:,3),'-r')
            hold on 
            plot(Unconstrained_data(:,1),Unconstrained_data(:,4),'-g')
            xlabel('Iteration #');ylabel('Number of Patches Developed')
            title({'Unconstrained Regulated Conventional';'Number of Patches Developed vs. Iteration #'})
            legend('Mussel','Finfish','Kelp','location','NorthWest')
            set(gcf,'color','white'); 
            axis tight 
            axis square
         % Constrained
            figure
            plot(Constrained_data(:,1),Constrained_data(:,2),'-b')
            hold on 
            plot(Constrained_data(:,1),Constrained_data(:,3),'-r')
            hold on 
            plot(Constrained_data(:,1),Constrained_data(:,4),'-g')
            xlabel('Iteration #');ylabel('Number of Patches Developed')
            title({'Constrained Regulated Conventional';'Number of Patches Developed vs. Iteration #'})
            legend('Mussel','Finfish','Kelp','location','NorthWest')
            set(gcf,'color','white'); 
            axis tight 
            axis square 
%% Pure Profit Driven Plans %%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% prompt_pure_profit='Run Pure Profit? = ';%input(prompt_pure_profit,'s');
str_pure_profit=Y;
if strcmp(str_pure_profit,'Y')
% Mussel impacted value 
    Mv_impacts_PP=M; [Y_M,I_M]=sort(Mv_impacts_PP,'descend');
% Finfish impacted value
    Fv_impacts_PP=F; [Y_F,I_F]=sort(Fv_impacts_PP,'descend');
% Kelp impacted value
    Kv_impacts_PP=K; [Y_K,I_K]=sort(Kv_impacts_PP,'descend');
% Use the cell FID as a common element between the farm types in the
% simulation 
    Aqua_FID=target_fid_fulldomain(V1+V2+V3>0);
    Aqua_FID_perm=Aqua_FID;
    
%% Pure Profit Unconstrained Model 
% Planning simulation loop 
    Suitability_Mat_PP=[Mv_impacts_PP Fv_impacts_PP Kv_impacts_PP];
    FID_Vector=Aqua_FID;
    FID_Vector_Reference=FID_Vector;
    Plan_initial=zeros(length(Aqua_dev_indices),1);
    Prev_plan=Plan_initial;
    tmp_plan=Prev_plan;
    itor=1;
    num_dev_patches=0;
    Unconstrained_Plan_Store_PP=[Plan_initial NaN(length(Aqua_dev_indices))];
    Unconstrained_data_PP=NaN(size(Unconstrained_Plan_Store_PP,2),4);
    while num_dev_patches<length(Aqua_dev_indices)
        itor=itor+1;
        [Y,I]=sort(Suitability_Mat_PP,1,'descend');
        [~,itor_aqua_type]=max(Y(1,:));
        Best_Cell_Y=Y(1,itor_aqua_type);
        Best_Cell_I=I(1,itor_aqua_type);
        Best_Cell_FID=FID_Vector(Best_Cell_I);
        [~,IA_1,~]=intersect(FID_Vector_Reference,Best_Cell_FID);% Identify the most suitable cell FID
        tmp_plan=Prev_plan;
        tmp_plan(IA_1)=itor_aqua_type;
        Unconstrained_Plan_Store_PP(:,itor)=tmp_plan;
        [~,IA_2,~]=intersect(FID_Vector,Best_Cell_FID);
        FID_Vector(IA_2)=[];
        Suitability_Mat_PP(IA_2,:)=[];
        num_dev_patches=sum(tmp_plan>0);
        disp(['iteration = ',num2str(itor)])
        disp(['Mussel = ',num2str(sum(tmp_plan==1)),' Finfish = ',num2str(sum(tmp_plan==2)),' Kelp = ',num2str(sum(tmp_plan==3))])
        if itor_aqua_type==1
         [~,IA_M_check,~]=intersect(target_fid_fulldomain(V1),FID_Vector_Reference(IA_1));
         if isempty(IA_M_check)==1
             disp('Error not developable cell for Mussels')
         end
        elseif itor_aqua_type==2
         [~,IA_F_check,~]=intersect(target_fid_fulldomain(V2),FID_Vector_Reference(IA_1));
         if isempty(IA_F_check)==1
             disp('Error not developable cell for Finfish')
         end  
        elseif itor_aqua_type==3
         [~,IA_K_check,~]=intersect(target_fid_fulldomain(V3),FID_Vector_Reference(IA_1));
         if isempty(IA_K_check)==1
             disp('Error not developable cell for Kelp')
         end
        end
        Unconstrained_data_PP(itor,1)=itor;Unconstrained_data_PP(itor,2)=sum(tmp_plan==1);Unconstrained_data_PP(itor,3)=sum(tmp_plan==2);
        Unconstrained_data_PP(itor,4)=sum(tmp_plan==3);
        Prev_plan=tmp_plan;
        disp(itor_aqua_type)
    end

        Unconstrained_Values_Static_PP=NaN(size(Unconstrained_Plan_Store_PP,2),7);
        for index=1:size(Unconstrained_Plan_Store_PP,2);
            X=zeros(6425,1);X(Aqua_dev_indices)=Unconstrained_Plan_Store_PP(:,index);
            % Mussel 
                Unconstrained_Values_Static_PP(index,1)=sum(Ms(Unconstrained_Plan_Store_PP(:,index)==1));
            % Finfish
                Unconstrained_Values_Static_PP(index,2)=sum(Fs(Unconstrained_Plan_Store_PP(:,index)==2));
            % Kelp
                Unconstrained_Values_Static_PP(index,3)=sum(Ks(Unconstrained_Plan_Store_PP(:,index)==3));
            % Halibut 
                Unconstrained_Values_Static_PP(index,4)=sum(Halibut_10yrNPVi(X==0))/sum(Halibut_10yrNPVi);
            % Viewshed
                Unconstrained_Values_Static_PP(index,5)=1-(sum((V_F.*(X==2))+(V_KM.*[(X==1)+(X==3)]))/sum(V));
            % Benthic 
                Unconstrained_Values_Static_PP(index,6)=sum(Bs(Unconstrained_Plan_Store_PP(:,index)~=2));
            % Disease 
                Unconstrained_Values_Static_PP(index,7)=sum(Ds(Unconstrained_Plan_Store_PP(:,index)~=2));
        end
%         load('Tuner_save_condensed_April_8_2016')
%         load('Tuner_save_Jun_22_2015_15_28','target_fid_fulldomain','Aqua_dev_indices','P1_full_scaled_obj','P2_full_scaled_obj','P3_full_scaled_obj','P4_aqua_patches_scaled_obj','Viewshed_impacts_finfish_scaled_obj',...
%          'Viewshed_impacts_kelpmussel_scaled_obj','P6_aqua_patches_scaled_obj','disease_min','disease_max','Fish_NPVi_indices','SQ_1fishable_0notfishable_for_each_soft_depth_patch_ORIGINAL'...
%          ,'Nij_initial_msy_TS','x0_PPUE_msy_TS','Fsum_msy','mussel_NPVi','Fish_NPVi','kelp_NPVi','P4_raw','Halibut_max','P4_full_dev','Viewshed_impacts_finfish','Viewshed_impacts_kelpmussel'...
%          ,'enviro_impact','target_fid_hab_depth','Wij','numpatches','max_age','T','delta','age_legal','age_mature','Dii','alphaCR','beta_i','habitat_area_i','theta','price','gamma',...
%         'distance_to_port_for_each_soft_depth_patch','Mii','habitat_area_i','age_move','Run_full_1_Run_dummy_0','dummy_indices','discount_rate_iy','disease_connect_matrix','x0_PPUE_initial','x0_PPUE_msy_TS')
        warning('off','all')
        gammatmp=gamma;
        Unconstrained_Dynamic_values_PP=NaN(size(Unconstrained_Plan_Store_PP,2),7);
        for indexEF=1:length(Unconstrained_Dynamic_values_PP)
            
        disp(['Unconstrained Iteration #',num2str(indexEF)])
        
        X=Unconstrained_Plan_Store_PP(:,indexEF);
        X_domain=zeros(size(target_fid_fulldomain));
        X_domain(Aqua_dev_indices)=X;
        [P1,P2,P3,P4,P5,P6,P7,P1i,P2i,P3i,P4i,P5i,P6i,P7i,mussel_NPVi_tmp,mussel_NPV_tmp,Fish_NPVi_tmp,Fish_NPV_tmp,kelp_NPVi_tmp,kelp_NPV_tmp,Halibut_NPVi_tmp,Halibut_NPV_tmp,Halibut_NPV_scaled_tmp,Halibut_Yi_tmp,halibut_dynamic_scaled_tmp,Viewshed_raw_valuei_fin_tmp,Viewshed_raw_valuei_kelpmussel_tmp,Viewshed_raw_value_fin_tmp,Viewshed_raw_value_kelpmussel_tmp,enviro_impacti_tmp,enviro_impact_tmp,disease_rawi,disease_raw_sum]=Sector_calcs_wrt_X(X,target_fid_fulldomain,Aqua_dev_indices,P1_full_scaled_obj,P2_full_scaled_obj,P3_full_scaled_obj,P4_aqua_patches_scaled_obj,Viewshed_impacts_finfish_scaled_obj,...
        Viewshed_impacts_kelpmussel_scaled_obj,P6_aqua_patches_scaled_obj,disease_min,disease_max,Fish_NPVi_indices,SQ_1fishable_0notfishable_for_each_soft_depth_patch_ORIGINAL,Nij_initial_msy_TS,x0_PPUE_msy_TS,Fsum_msy,mussel_NPVi,Fish_NPVi,kelp_NPVi,P4_raw,Halibut_max,P4_full_dev,Viewshed_impacts_finfish,Viewshed_impacts_kelpmussel,enviro_impact,target_fid_hab_depth,Wij,numpatches,max_age,T,delta,age_legal,age_mature,Dii,alphaCR,beta_i,habitat_area_i,theta,price,gammatmp,distance_to_port_for_each_soft_depth_patch,Mii,age_move,Run_full_1_Run_dummy_0,dummy_indices,discount_rate_iy,disease_connect_matrix);
        
        Unconstrained_Dynamic_values_PP(indexEF,1)=P1; 
        Unconstrained_Dynamic_values_PP(indexEF,2)=P2;
        Unconstrained_Dynamic_values_PP(indexEF,3)=P3;
        Unconstrained_Dynamic_values_PP(indexEF,4)=Halibut_Yi_tmp./Halibut_max;disp(num2str(Halibut_Yi_tmp./Halibut_max))
        Unconstrained_Dynamic_values_PP(indexEF,5)=1-(sum((V_F(Aqua_dev_indices).*(X==2))+(V_KM(Aqua_dev_indices).*[(X==1)+(X==3)]))/sum(V));
        Unconstrained_Dynamic_values_PP(indexEF,6)=P6;
        Unconstrained_Dynamic_values_PP(indexEF,7)=P7; 
        
        end
        save('Pure_Profit_Unconstrained_Dynamic_Values','Unconstrained_Dynamic_values_PP')
%% Pure Profit Constrained Model 
% Load floating variables 
    Suitability_Mat_PP=[Mv_impacts_PP Fv_impacts_PP Kv_impacts_PP];
    FID_Vector=Aqua_FID;
    FID_Vector_Reference=FID_Vector;
    Plan_initial=zeros(length(Aqua_dev_indices),1);
    Prev_plan=Plan_initial;
    tmp_plan=Prev_plan;
    itor=1;
    num_dev_patches=0;
    Constrained_Plan_Store_PP=NaN(length(Aqua_dev_indices));
    Constrained_data_PP=NaN(size(Constrained_Plan_Store_PP,2),4);
    while num_dev_patches<length(Aqua_dev_indices)
        n_aqua_types=find(sum(Suitability_Mat_PP>0,1)>0);
        for index=1:length(n_aqua_types)
            itor_aqua_type=n_aqua_types(index);
            itor=itor+1;
            [Y,I]=sort(Suitability_Mat_PP,'descend');
            Best_Cell_Y=Y(1,itor_aqua_type);
            Best_Cell_I=I(1,itor_aqua_type);
            Best_Cell_FID=FID_Vector(Best_Cell_I);
            [~,IA_1,~]=intersect(FID_Vector_Reference,Best_Cell_FID);% Identify the most suitable cell FID
            tmp_plan=Prev_plan;
            tmp_plan(IA_1)=itor_aqua_type;
            Constrained_Plan_Store_PP(:,itor)=tmp_plan;
            [~,IA_2,~]=intersect(FID_Vector,Best_Cell_FID);
            FID_Vector(IA_2)=[];
            Suitability_Mat_PP(IA_2,:)=[];
            num_dev_patches=sum(tmp_plan>0);
            disp(['iteration = ',num2str(itor)])
            disp(['Mussel = ',num2str(sum(tmp_plan==1)),' Finfish = ',num2str(sum(tmp_plan==2)),' Kelp = ',num2str(sum(tmp_plan==3))])
            if itor_aqua_type==1
             [~,IA_M_check,~]=intersect(target_fid_fulldomain(V1),FID_Vector_Reference(IA_1));
             if isempty(IA_M_check)==1
                 disp('Error not developable cell for Mussels')
             end
            elseif itor_aqua_type==2
             [~,IA_F_check,~]=intersect(target_fid_fulldomain(V2),FID_Vector_Reference(IA_1));
             if isempty(IA_F_check)==1
                 disp('Error not developable cell for Finfish')
             end  
            elseif itor_aqua_type==3
             [~,IA_K_check,~]=intersect(target_fid_fulldomain(V3),FID_Vector_Reference(IA_1));
             if isempty(IA_K_check)==1
                 disp('Error not developable cell for Kelp')
             end
            end
            Constrained_data_PP(itor,1)=itor;Constrained_data_PP(itor,2)=sum(tmp_plan==1);Constrained_data_PP(itor,3)=sum(tmp_plan==2);
            Constrained_data_PP(itor,4)=sum(tmp_plan==3);
            Prev_plan=tmp_plan;
            disp(itor_aqua_type)
        end
    end
    Constrained_Plan_Store_PP(:,1)=zeros(length(Aqua_dev_indices),1);
    Constrained_Values_Static_PP=NaN(size(Constrained_Plan_Store_PP,2),7);
    for index=1:size(Constrained_Plan_Store_PP,2);
        X=zeros(6425,1);X(Aqua_dev_indices)=Constrained_Plan_Store_PP(:,index);
        % Mussel 
            Constrained_Values_Static_PP(index,1)=sum(Ms(Constrained_Plan_Store_PP(:,index)==1));
        % Finfish
            Constrained_Values_Static_PP(index,2)=sum(Fs(Constrained_Plan_Store_PP(:,index)==2));
        % Kelp
            Constrained_Values_Static_PP(index,3)=sum(Ks(Constrained_Plan_Store_PP(:,index)==3));
        % Halibut 
            Constrained_Values_Static_PP(index,4)=sum(Halibut_10yrNPVi(X==0))/sum(Halibut_10yrNPVi);
        % Viewshed
            Constrained_Values_Static_PP(index,5)=1-(sum((V_F.*(X==2))+(V_KM.*[(X==1)+(X==3)]))/sum(V));
        % Benthic 
            Constrained_Values_Static_PP(index,6)=sum(Bs(Constrained_Plan_Store_PP(:,index)~=2));
        % Disease 
            Constrained_Values_Static_PP(index,7)=sum(Ds(Constrained_Plan_Store_PP(:,index)~=2));
    end
    %% Run Dynamic Model 
%         load('Tuner_save_condensed_April_8_2016')
%         load('Tuner_save_Jun_22_2015_15_28','target_fid_fulldomain','Aqua_dev_indices','P1_full_scaled_obj','P2_full_scaled_obj','P3_full_scaled_obj','P4_aqua_patches_scaled_obj','Viewshed_impacts_finfish_scaled_obj',...
%         'Viewshed_impacts_kelpmussel_scaled_obj','P6_aqua_patches_scaled_obj','disease_min','disease_max','Fish_NPVi_indices','SQ_1fishable_0notfishable_for_each_soft_depth_patch_ORIGINAL'...
%         ,'Nij_initial_msy_TS','x0_PPUE_msy_TS','Fsum_msy','mussel_NPVi','Fish_NPVi','kelp_NPVi','P4_raw','Halibut_max','P4_full_dev','Viewshed_impacts_finfish','Viewshed_impacts_kelpmussel'...
%         ,'enviro_impact','target_fid_hab_depth','Wij','numpatches','max_age','T','delta','age_legal','age_mature','Dii','alphaCR','beta_i','habitat_area_i','theta','price','gamma',...
%         'distance_to_port_for_each_soft_depth_patch','Mii','habitat_area_i','age_move','Run_full_1_Run_dummy_0','dummy_indices','discount_rate_iy','disease_connect_matrix','x0_PPUE_initial','x0_PPUE_msy_TS')
        warning('off','all')
        gammatmp=gamma;
        Constrained_Dynamic_values_PP=NaN(size(Constrained_Plan_Store_PP,2),7);
        X_iter=Constrained_Plan_Store_PP;
        for indexEF=1:length(Constrained_Dynamic_values_PP)
        disp(['Constrained Iteration #',num2str(indexEF)])
        X=X_iter(:,indexEF);
        X_domain=zeros(size(target_fid_fulldomain));
        X_domain(Aqua_dev_indices)=X;
        [P1,P2,P3,P4,P5,P6,P7,P1i,P2i,P3i,P4i,P5i,P6i,P7i,mussel_NPVi_tmp,mussel_NPV_tmp,Fish_NPVi_tmp,Fish_NPV_tmp,kelp_NPVi_tmp,kelp_NPV_tmp,Halibut_NPVi_tmp,Halibut_NPV_tmp,Halibut_NPV_scaled_tmp,Halibut_Yi_tmp,halibut_dynamic_scaled_tmp,Viewshed_raw_valuei_fin_tmp,Viewshed_raw_valuei_kelpmussel_tmp,Viewshed_raw_value_fin_tmp,Viewshed_raw_value_kelpmussel_tmp,enviro_impacti_tmp,enviro_impact_tmp,disease_rawi,disease_raw_sum]=Sector_calcs_wrt_X(X,target_fid_fulldomain,Aqua_dev_indices,P1_full_scaled_obj,P2_full_scaled_obj,P3_full_scaled_obj,P4_aqua_patches_scaled_obj,Viewshed_impacts_finfish_scaled_obj,...
        Viewshed_impacts_kelpmussel_scaled_obj,P6_aqua_patches_scaled_obj,disease_min,disease_max,Fish_NPVi_indices,SQ_1fishable_0notfishable_for_each_soft_depth_patch_ORIGINAL,Nij_initial_msy_TS,x0_PPUE_msy_TS,Fsum_msy,mussel_NPVi,Fish_NPVi,kelp_NPVi,P4_raw,Halibut_max,P4_full_dev,Viewshed_impacts_finfish,Viewshed_impacts_kelpmussel,enviro_impact,target_fid_hab_depth,Wij,numpatches,max_age,T,delta,age_legal,age_mature,Dii,alphaCR,beta_i,habitat_area_i,theta,price,gammatmp,distance_to_port_for_each_soft_depth_patch,Mii,age_move,Run_full_1_Run_dummy_0,dummy_indices,discount_rate_iy,disease_connect_matrix);
        
        Constrained_Dynamic_values_PP(indexEF,1)=P1; 
        Constrained_Dynamic_values_PP(indexEF,2)=P2;
        Constrained_Dynamic_values_PP(indexEF,3)=P3;
        Constrained_Dynamic_values_PP(indexEF,4)=Halibut_Yi_tmp./Halibut_max;disp(num2str(Halibut_Yi_tmp./Halibut_max))
        Constrained_Dynamic_values_PP(indexEF,5)=1-(sum((V_F(Aqua_dev_indices).*(X==2))+(V_KM(Aqua_dev_indices).*[(X==1)+(X==3)]))/sum(V));
        Constrained_Dynamic_values_PP(indexEF,6)=P6;
        Constrained_Dynamic_values_PP(indexEF,7)=P7; 
        
        end
        save('Pure_Profit_Constrained_Dynamic_Values','Constrained_Dynamic_values_PP')
    %% Save variables
        filename = 'C:\Users\Joel\Desktop\Thesis YTD\Code\MSP Planning Results April 2016';
        cd(filename)
        save('Pure_Profit_Conventional_Results_April','Unconstrained_Plan_Store','Unconstrained_Values_Static','Constrained_Plan_Store','Constrained_Values_Static')
        
        csvwrite('Pure_Profit_Unconstrained_Plans_April.csv',Unconstrained_Plan_Store_PP)
        csvwrite('Pure_Profit_Unconstrained_Static_Values_April.csv',Unconstrained_Values_Static_PP)
        
        csvwrite('Pure_Profit_Constrained_Plans_April.csv',Constrained_Plan_Store_PP)
        csvwrite('Pure_Profit_Constrained_Static_Values_April.csv',Constrained_Values_Static_PP)
        
        csvwrite('Pure_Profit_Unconstrained_Dynamic_Values_April.csv',Unconstrained_Dynamic_values_PP)
        csvwrite('Pure_Profit_Constrained_Dynamic_Values_April.csv',Constrained_Dynamic_values_PP)
        
        filename = 'C:\Users\Joel\Desktop\Thesis YTD\Code';
        cd(filename)
    %% Plot conventional management figures 
        % Unconstrained 
            figure
            plot(Unconstrained_data_PP(:,1),Unconstrained_data_PP(:,2),'-b')
            hold on 
            plot(Unconstrained_data_PP(:,1),Unconstrained_data_PP(:,3),'-r')
            hold on 
            plot(Unconstrained_data_PP(:,1),Unconstrained_data_PP(:,4),'-g')
            xlabel('Iteration #');ylabel('Number of Patches Developed')
            title({'Pure Profit Unconstrained Regulated Conventional';'Number of Patches Developed vs. Iteration #'})
            legend('Mussel','Finfish','Kelp','location','NorthWest')
            set(gcf,'color','white'); 
            axis tight 
            axis square
         % Constrained
            figure
            plot(Constrained_data_PP(:,1),Constrained_data_PP(:,2),'-b')
            hold on 
            plot(Constrained_data_PP(:,1),Constrained_data_PP(:,3),'-r')
            hold on 
            plot(Constrained_data_PP(:,1),Constrained_data_PP(:,4),'-g')
            xlabel('Iteration #');ylabel('Number of Patches Developed')
            title({'Pure Profit Constrained Regulated Conventional';'Number of Patches Developed vs. Iteration #'})
            legend('Mussel','Finfish','Kelp','location','NorthWest')
            set(gcf,'color','white'); 
            axis tight 
            axis square 
end