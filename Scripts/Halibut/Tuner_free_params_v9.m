close all
clear all
disp('------------------------')
%Nomenclature of subscripts to parameters used below
% i = patch = rows
% j = age = columns
tic;
%% Settings
warning('off'); %no warning e.g. when divide by zero
%Set defaults for making the figures prettyse
set(0,'DefaultTextFontSize',35)
set(0,'DefaultAxesFontSize',35)
set(0,'DefaultLineLineWidth',6)
set(0,'DefaultAxesLineWidth',2)
set(0,'DefaultSurfaceLineWidth',2)
%% globals
pathtool
w = waitforbuttonpress;
if w == 0
    disp('Button click')
else
    disp('Key press')
end
load 'file_dir_params'
global alphaCR habitat_area_i numpatches max_age Bij age_mature Dii K age_legal Ei delta Nij Mii age_move Wij TargetBiomass beta_i...
    Rmax CRgoal Nij_initial theta price Theta_Density_prop K_legal_virgin_abundance PostHstockdensity_i PostHstockdensity_kvirgin...
    Fsum Bi_legal Kvirgin_abundance SSB_msy_target Fsum_msy_guess gamma distance_to_port_for_each_soft_depth_patch unique_block_id...
    blockid_for_each_soft_depth_patch PacCoastFisheryGIS_kglandings_by_block_Com_plus_Rec_sum2one Nij_initial_virgin...
    actual_comm_plus_rec_kglandings CRtarget Fsum_UB SQ_1fishable_0notfishable_for_each_soft_depth_patch SQ_1trawling_0notrawling_for_each_soft_detph_patch...
    qi depth_i habitat_area_i_raw depth_bins runsetparams_1full col_soft_hab_area study_area_polygons_PacCoastFisheryGIS_NUM study_area_polygons_PacCoastFisheryGIS_TXT_noheader...
    col_row_labels cells_with_soft_bottom_NUM connect_matrix Halibut_adult_movement_v3_1full halibut_max_depth max_possible_depth Fsum_opt...
    habitat_area_i_fulldomain depth_i_fulldomain Nij_Fsumchoice nearestdist_to_outfall omega study_area_polygons_PacCoastFisheryGIS_TXT_noheader_soft_depthI...
    cells_with_hab_NUM tol_simulation T x0_PPUE_initial discount_rate
%% Code task choices - CHOOSE
load0_do1_tuning=1; %Tune the free parameters, or load pre-tuned parameters
%The below options only work if load0_do1_tuning=0
Plot_params1_0N_1Y=1; %Plot the spatial and demographic parameters? Note: plotting the adult movement and larval dispersal kernels is slow
Verify_model_0N_1Y=1; %Run core and fleet models under virgin and MSY conditions and compare/critique results
runtimeseries_0N_1Y=1; %Run model over finite time horizon and record spatially and temporally explicit results (e.g., for NPV calculation)
% run_dyn_eval_MSPsolns_0N_1Y=1; %Evaluate dynamic response to each MSP aquaculture solution

%TS model equilibrium params - DO NOT CHANGE
T=500; %number of years to run
discount_rate=0; %discount rate

%% Tune unnknown (free) parameters
    toc_tuner=toc;
if load0_do1_tuning==1
    disp('Tune free parameters')

    % Set depth
    max_possible_depth=-90; %Haugen 1990 Figure 2
    halibut_max_depth=max_possible_depth;

    runsetparams_1full=1; %run full code, including data loading
    Set_params3
    runsetparams_1full=0; %don't need full code to run again

    %initial tuning guesses
    Rmax_guess=0.00068751;
    alphaCR_guess = 24.4282;
    Fsum_msy_guess=2.043675079230346e+03;
    gamma_guess=6.0708e-05;

    %Psi reduces habitat quality wrt depth (Haugen 1990)
    psi_Haugen=-0.233; psi=psi_Haugen;
    %given psi, find depth where hab prop drops to zero. Set max depth
    %at just above that level (so as to avoid modeling zero habitat patches)
    depth_bins=linspace(0,max_possible_depth,100)';
    habitat_area_propreduct_wrt_depth=habitat_area_i_wrt_habitat_prop_reduction_wrt_depth(psi,depth_bins,ones(size(depth_bins)));
    tmp1=find(habitat_area_propreduct_wrt_depth==0);
        if min(habitat_area_propreduct_wrt_depth)>0
            halibut_max_depth=max_possible_depth;
        else
            halibut_max_depth=ceil(depth_bins(tmp1(1)-1));
        end
    %load params, some of which are wrt halibut_max_depth:
    Set_params3
    %now set habitat area wrt to the chosen patches and psi
    habitat_area_i=habitat_area_i_wrt_habitat_prop_reduction_wrt_depth(psi,depth_i,habitat_area_i_raw);

    %calculate adult movement
    Halibut_adult_movement_v4
    Mii=mij; %adult dispersal matrix

    %initial conditions
    Nij_initial_virgin=ones(numpatches,max_age); %matrix of initial number of individuals age j (col) in each patch i (row); placeholder
    Nij_initial=Nij_initial_virgin; %Set initial conditions for model. Note: setting to virgin means the model will not run fast (because has to harvest down from virgin every time), but it seems to help with the tuning
    tol_tuner=0.01; %proportion difference allowed make when tuning params
    toc_fulltuning=toc;

    %starting free param values (placeholders)
    alphaCR=alphaCR_guess;
    gamma=gamma_guess; %placeholder

    diff_tuner=1; %placeholder
    Ti=0;
%     break
    while diff_tuner>tol_tuner
        Ti=Ti+1;
        disp(['Tuning iteration: ',num2str(Ti)])

        %tune Rmax so that at Fsum=Fsum_msy, sum(Yi)=Actual yield
        toc_Rmax=toc;
        TolFun_value=1e-4*numpatches;
        options=optimset('Display','off','TolFun',TolFun_value);
        [X,FVAL,EXITFLAG]= fmincon(@tuneRmax_with_Fsumopt,Rmax_guess,[],[],[],[],1e-5,0.01,[],options);
        diff_Rmax=abs(X-Rmax_guess)/X;
        disp(['Tuned Rmax = ',num2str(X),'; Took ',num2str((toc-toc_Rmax)/60),' minutes. Flag = ',num2str(EXITFLAG),])
        Rmax_guess=X;
        Rmax=X;

        %tune alphaCR so that CR=CRtarget under virgin conditions
        toc_alphaCR=toc;
        TolX_value=1e-4*numpatches; %Set tolerance for alphaCR search. 1e-4 is the default TolX.
        TolFun_value=TolX_value;
        options=optimset('Display','off','TolX',TolX_value,'TolFun',TolFun_value);
        [X,FVAL,EXITFLAG]= fmincon(@tuneAlphaCRtarget,alphaCR_guess,[],[],[],[],5,50,[],options);
        disp(['Tuned alphaCR = ',num2str(X),'; Took ',num2str((toc-toc_alphaCR)/60),' minutes. Flag = ',num2str(EXITFLAG),])
        diff_alphaCR=abs(X-alphaCR_guess)/X; %calc diff from last iteration
        alphaCR_guess=X; %set alphaCR guess for next iteration (if it happens)
        alphaCR=X; %set alphaCR for below evaluation
        %Recalculate beta and theta givenn Rmax and alphaCR choices:
        beta_i=alphaCR./(Rmax.*habitat_area_i); % calc beta wrt alphaCR
        %run w/o harvest so can calculate theta
        Ei=zeros(numpatches,1); %set no harvest to can determine PostHstockdensity_kvirgin for calculating theta
        Nij_initial=Nij_initial_virgin; %set initial conditions (note: Nij_initial_virgin is reset in Rmax and psi tuners)
        theta=0; %placeholder
        BasicSpatialHarvestModel_COREv1_TS %run model to virgin K
        PostHstockdensity_kvirgin_i=PostHstockdensity_i; %stock biomass density per patch
        PostHstockdensity_kvirgin=(sum(PostHstockdensity_kvirgin_i.*habitat_area_i))/sum(habitat_area_i); %average across all habitat patches of stock biomass density per unit area habitat
        theta=PostHstockdensity_kvirgin*Theta_Density_prop*price; %set theta

        %re-find Fsum_msy at alphaCR=X;
        Fsum_UB=sum(habitat_area_i)*4.5e-6; %set max effort available for search
        TolX_value=1e-4*numpatches; %Set tolerance for Fsum search. 1e-4 is the default TolX.
        TolFun_value=TolX_value;
        % We have so many patches and thus Fsum is so large (and fracional changes in Fsum don't matter much), so make Fsum tolerance large
        options=optimset('Display','off','TolX',TolX_value,'TolFun',TolFun_value);
        [Fsum_opt negprofit flag]=fmincon(@OptFsum,Fsum_msy_guess,[],[],[],[],0,Fsum_UB,[],options);
        if Fsum_opt==Fsum_UB; disp(['Fsum_opt being constrainted by Fsum_UB=',num2str(Fsum_UB),' !!!%*^%$##@']); end
        Fsum_msy_guess=Fsum_opt;

        %Now run with Fsum=Fsum_msy
        Fsum=Fsum_opt;
        Nij_initial=Nij_Fsumchoice; %based on Fsum opt from Rmax tuning
        BasicSpatialHarvestModel_COREv1_fleet_TS

        %Tune distance-from-port cost param gamma so that spatial distribution of
        %landings at Fsum_msy matches the empirical CDFW block data
        toc_gamma=toc; %gamma optimization start time
        gamma_LB=0; %lower possible bound for gamma
        distance_max_amplifier=100; % (1+distance_max_amplifier)=proportional increase in cost to the farthest patch (i.e., 2 =twice as expensive as patch at port)
        gamma_UB=(distance_max_amplifier./max(distance_to_port_for_each_soft_depth_patch)); %gamma upper bound
        options=optimset('Display','off');
        [X,FVAL,EXITFLAG]=fmincon(@tunegamma_wrt_spatialEi,gamma_guess,[],[],[],[],gamma_LB,gamma_UB,[],options);
        diff_gamma=abs(X-gamma_guess)/X;
        if gamma_guess==0 && X==0; diff_gamma=0; end %then no change in this tuning iteration
        gamma_guess=X;
        disp(['Tuned gamma = ',num2str(X),'; SS diffs in distribution = ',num2str(FVAL),'; Took ',num2str((toc-toc_gamma)/60),' minutes. Flag = ',num2str(EXITFLAG),])
        gamma=X;

        %When this gets small enough tuner while....end loop is complete and tuning is complete
        diff_tuner=max([diff_Rmax diff_alphaCR diff_gamma]);
    end

    Fsum_msy=Fsum_opt; %Fsum; %store!
    %Run fleet model with tuned params. First do no harvest to get virgin outputs:
    Ei=zeros(numpatches,1); %set no harvest to can determine PostHstockdensity_kvirgin for calculating theta
    theta=0; %placeholder
    BasicSpatialHarvestModel_COREv1_TS %run model to virgin K
    CRvirgin=CR; %get CR
    PostHstockdensity_kvirgin_i=PostHstockdensity_i; %stock biomass density per patch
    PostHstockdensity_kvirgin=(sum(PostHstockdensity_kvirgin_i.*habitat_area_i))/sum(habitat_area_i); %average across all habitat patches of stock biomass density per unit area habitat
    theta=PostHstockdensity_kvirgin*Theta_Density_prop*price; %set theta
    SSBvirgin=SSB_postharvest; %get virgin SSB
    Nij_initial_virgin=Nij; %get Nij virgin
    %Now run model at Fsum=Fsum_msy
    Fsum=Fsum_msy;
    BasicSpatialHarvestModel_COREv1_fleet_TS
    SSB_wrt_Fsum_msy=SSB_postharvest;
    Yield_Fsum_msy=sum(Yi);
    Ei_wrt_Fsum_msy=Ei;
    Nij_initial_msy_TS=Nij; %Use this input for MSP analysis
    x0_PPUE_msy_TS=x0_PPUE; %Use this input for MSP analysis
    %Display tuned values:
    disp(['Final tuned values: Rmax = ',num2str(Rmax),'; alphaCR = ',num2str(alphaCR),'; gamma = ',num2str(gamma),'; psi = ',num2str(psi)...
        '; CRvirgin = ',num2str(CRvirgin),'; Fsum_msy = ',num2str(Fsum_msy),'; Ei_wrt_Fsum_msy = ',num2str(min(Ei_wrt_Fsum_msy)),'-',num2str(max(Ei_wrt_Fsum_msy)),...
        '; delta=',num2str(delta(1)),'; Took ',num2str(Ti),' iterations and ',num2str((toc-toc_tuner)/60/60),' hours; diff_tuner=',num2str(diff_tuner)])
    disp(['Tuned outputs: actual_comm_plus_rec_kglandings/Yield_Fsum_msy = ',num2str(actual_comm_plus_rec_kglandings/Yield_Fsum_msy),'; SSB_wrt_Fsum_msy/SSBvirgin=',num2str(SSB_wrt_Fsum_msy/SSBvirgin),])
    %clear toggles so not saved so don't overwrite toggle choice when loading tuned_params
    clear load0_do1_tuning Plot_params1_0N_1Y Verify_model_0N_1Y runtimeseries_0N_1Y
    %save all the tuned parameter results
    tuned_params=Fsum;
    save tuned_params
else %do not tune params...just load them from previous tuning
    load('tuned_params')
    disp(['Loaded tuned values: Rmax = ',num2str(Rmax),'; alphaCR = ',num2str(alphaCR),'; gamma = ',num2str(gamma),'; psi = ',num2str(psi)])
    discount_factor=1./((1+discount_rate).^[1:T]);
    discount_rate_iy=repmat(discount_factor,numpatches,1);
end
% break
%% Verify and plot results based on tuned params:
if Plot_params1_0N_1Y==1
    Plot_params1
end
if Verify_model_0N_1Y==1
%verify that loaded params are generating properly tuned model results -
%Equilibrium simulations:
% I. set harvest to zero and run virgin equil results
    %First do with core model that uses raw spatial fishing effort levels (Ei)
    Ei=zeros(size(Ei_wrt_Fsum_msy));
    Nij_initial=Nij_initial_virgin;
    BasicSpatialHarvestModel_COREv1_TS
    SSBvirgin=SSB_postharvest; %get virgin SSB
    Yield_virgin=sum(Yi); %yield
    disp(['Virgin equilibrium model outputs (Ei model): CR = ',num2str(CR),'; SSBvirgin=',num2str(SSBvirgin),' kg; Yield=',num2str(Yield_virgin)])
    panelplot_biomass_yield_equilprofit_npv; suptitle('Calibrated model output: Virgin [Core model]')
    %Second do with fleet model that uses Fsum (results should be the same!)
    Fsum=0;
    BasicSpatialHarvestModel_COREv1_fleet_TS
    SSBvirgin=SSB_postharvest; %get virgin SSB
    SSBivirgin=sum(SSBij_postharvest,2);
    Yield_virgin=sum(Yi); %yield
    disp(['Virgin equilibrium model outputs (Fsum model): CR = ',num2str(CR),'; SSBvirgin=',num2str(SSBvirgin),' kg; Yield=',num2str(Yield_virgin)])
    panelplot_biomass_yield_equilprofit_npv; suptitle('Calibrated model output: Virgin [Fleet model]')
% II. run model at spatial harvest levels that produce MSY
    %First do with core model that uses raw spatial fishing effort levels (Ei)
    Ei=Ei_wrt_Fsum_msy;
    Nij_initial=Nij_initial_msy_TS;
    x0_PPUE_initial=x0_PPUE_msy_TS;
    BasicSpatialHarvestModel_COREv1_TS
    SSB_wrt_Fsum_msy=SSB_postharvest;
    Yield_Fsum_msy=sum(Yi);
    disp(['MSY equilibrium model outputs (Ei model): Ei_wrt_Fsum_msy = ',num2str(min(Ei_wrt_Fsum_msy)),'-',num2str(max(Ei_wrt_Fsum_msy)),'; sum(Ei_wrt_Fsum_msy)/Fsum_msy=',num2str(sum(Ei_wrt_Fsum_msy)/Fsum_msy)])
    disp(['   actual_comm_plus_rec_kglandings/Yield_Fsum_msy = ',num2str(actual_comm_plus_rec_kglandings/Yield_Fsum_msy),'; SSB_wrt_Fsum_msy/SSBvirgin=',num2str(SSB_wrt_Fsum_msy/SSBvirgin)])
    panelplot_biomass_yield_equilprofit_npv; suptitle('Calibrated model output: MSY [Core model]')
    %the the observed and simulated landings in each CDFW fishing block:
    plot_CDFWblock_observed_model_data
    %Second do with fleet model that uses Fsum (results should be the same!)
    Fsum=Fsum_msy;
    BasicSpatialHarvestModel_COREv1_fleet_TS
    SSB_wrt_Fsum_msy=SSB_postharvest;
    SSBi_wrt_Fsum_msy=sum(SSBij_postharvest,2);
    Yield_Fsum_msy=sum(Yi);
    Yi_Fsum_msy=Yi;
    Payoffi_Fsum_msy=Payoff;
    disp(['MSY equilibrium model outputs (Fsum model): Fsum_msy=',num2str(Fsum_msy),'; Ei_wrt_Fsum_msy = ',num2str(min(Ei_wrt_Fsum_msy)),'-',num2str(max(Ei_wrt_Fsum_msy)),'; sum(Ei_wrt_Fsum_msy)/Fsum_msy=',num2str(sum(Ei_wrt_Fsum_msy)/Fsum_msy)])
    disp(['   actual_comm_plus_rec_kglandings/Yield_Fsum_msy = ',num2str(actual_comm_plus_rec_kglandings/Yield_Fsum_msy),'; SSB_wrt_Fsum_msy/SSBvirgin=',num2str(SSB_wrt_Fsum_msy/SSBvirgin)])
    panelplot_biomass_yield_equilprofit_npv; suptitle('Calibrated model output: MSY [Fleet model]')
    plot_map_results
    plot_TS_results
end
%% Near-term simulations
% close all
if runtimeseries_0N_1Y==1
    disp('Run model over finite time horizon and record spatially and temporally explicit results')
    %Choose below params
    T=10; %number of years to run (T=10)
    discount_rate=0.05; %discount rate (d=0.05)
    discount_factor=1./((1+discount_rate).^[1:T]);
    discount_rate_iy=repmat(discount_factor,numpatches,1);
    Nij_initial=Nij_initial_msy_TS; %initial stock status (patch and age class specific)
    x0_PPUE_initial=x0_PPUE_msy_TS;

    %Run with status quo MSY management and no Aquaculture
    disp('TS model with SQ MSY conditions...')
    Fsum=Fsum_msy; %TAE (fleet model)
    BasicSpatialHarvestModel_COREv1_fleet_TS
    plot_TS_results %Uncomment for plots
    %calc patch specific NPV of fishery
    tmp1=Yiy.*discount_rate_iy; %multiply year and patch specific values by discounted value
    Yi_NPV_MSY=sum(tmp1,2); %sum years together to calc NPV for each patch
    %create and save a vector of NPVs across entire domain of patches (not just halibut patches)
    Yi_fulldomain_NPV=zeros(size(filter_habitat_rows)); %vector of zeros for ALL patches in domain (not just halibut patches)
    Yi_fulldomain_NPV(filter_habitat_rows==1)=Yi_NPV_MSY; %fill in halibut patches with fishery values
    target_fid_fulldomain=study_area_polygons_PacCoastFisheryGIS_NUM(:,col_target_fid);
    xlswrite(strcat(output_data_dir,'Target_FID_and_Yi_fulldomain_NPV_at_MSY_noAqua.xlsx'),[target_fid_fulldomain Yi_fulldomain_NPV]);
end
exit;
