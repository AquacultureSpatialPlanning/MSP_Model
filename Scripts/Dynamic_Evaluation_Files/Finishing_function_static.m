    %% Finishing Function 
    function finish_out=Finishing_function_static(tasks_in)
    load('Tuner_save_Jun_22_2015_15_28','target_fid_fulldomain','Aqua_dev_indices','P1_full_scaled_obj','P2_full_scaled_obj','P3_full_scaled_obj','P4_aqua_patches_scaled_obj','Viewshed_impacts_finfish_scaled_obj',...
             'Viewshed_impacts_kelpmussel_scaled_obj','P6_aqua_patches_scaled_obj','disease_min','disease_max','Fish_NPVi_indices','SQ_1fishable_0notfishable_for_each_soft_depth_patch_ORIGINAL'...
             ,'Nij_initial_msy_TS','x0_PPUE_msy_TS','Fsum_msy','mussel_NPVi','Fish_NPVi','kelp_NPVi','P4_raw','Halibut_max','P4_full_dev','Viewshed_impacts_finfish','Viewshed_impacts_kelpmussel'...
             ,'enviro_impact','target_fid_hab_depth','Wij','numpatches','max_age','T','delta','age_legal','age_mature','Dii','alphaCR','beta_i','habitat_area_i','theta','price','gamma',...
            'distance_to_port_for_each_soft_depth_patch','Mii','habitat_area_i','age_move','Run_full_1_Run_dummy_0','dummy_indices','discount_rate_iy','disease_connect_matrix','x0_PPUE_initial','x0_PPUE_msy_TS')
    warning('off','all')
    gammatmp=gamma;
    % Load all data for the task specific structure 
        X_iter=tasks_in.Task_solutions;alpha_matrix_loop=tasks_in.Alphas;
        indexEF_tasks=tasks_in.Index_mat;
    Master_matrix_static=NaN(length(indexEF_tasks),29);
    for indexEF_iterations=1:length(tasks_in.Iteration_index)
            X=X_iter(:,indexEF_iterations);
            indexEF=indexEF_tasks(indexEF_iterations);
            [alpha1,alpha2,alpha3,alpha4,alpha5,alpha6,alpha7]=Tradeoff_solution_alpha_matrix(indexEF,alpha_matrix_loop); 
%             cd('C:\Users\jsteve13\Desktop\Static MSP\Static Solutions')
% %         if exist(['data_save_MSP_static_values_A1_', num2str(100*alpha1),'_A2_', num2str(100*alpha2),'_A3_', num2str(100*alpha3),'_A4_',num2str(100*alpha4),'_A5_', num2str(100*alpha5),'_A6_', num2str(100*alpha6),'_A7_', num2str(100*alpha7),'.mat'],'file')>0
% %             cd('C:\Users\jsteve13\Desktop\Static MSP')
% %             continue
% %         else 
%             cd('F:\Static MSP')
            [P1,P2,P3,P4,P5,P6,P7,P1i,P2i,P3i,P4i,P5i,P6i,P7i,mussel_NPVi_tmp,mussel_NPV_tmp,Fish_NPVi_tmp,Fish_NPV_tmp,kelp_NPVi_tmp,kelp_NPV_tmp,Halibut_NPVi_tmp,Halibut_NPV_tmp,Halibut_NPV_scaled_tmp,Halibut_Yi_tmp,halibut_dynamic_scaled_tmp,Viewshed_raw_valuei_fin_tmp,Viewshed_raw_valuei_kelpmussel_tmp,Viewshed_raw_value_fin_tmp,Viewshed_raw_value_kelpmussel_tmp,enviro_impacti_tmp,enviro_impact_tmp,disease_rawi,disease_raw_sum]=Sector_calcs_wrt_X(X,target_fid_fulldomain,Aqua_dev_indices,P1_full_scaled_obj,P2_full_scaled_obj,P3_full_scaled_obj,P4_aqua_patches_scaled_obj,Viewshed_impacts_finfish_scaled_obj,...
            Viewshed_impacts_kelpmussel_scaled_obj,P6_aqua_patches_scaled_obj,disease_min,disease_max,Fish_NPVi_indices,SQ_1fishable_0notfishable_for_each_soft_depth_patch_ORIGINAL,Nij_initial_msy_TS,x0_PPUE_msy_TS,Fsum_msy,mussel_NPVi,Fish_NPVi,kelp_NPVi,P4_raw,Halibut_max,P4_full_dev,Viewshed_impacts_finfish,Viewshed_impacts_kelpmussel,enviro_impact,target_fid_hab_depth,Wij,numpatches,max_age,T,delta,age_legal,age_mature,Dii,alphaCR,beta_i,habitat_area_i,theta,price,gammatmp,distance_to_port_for_each_soft_depth_patch,Mii,age_move,Run_full_1_Run_dummy_0,dummy_indices,discount_rate_iy,disease_connect_matrix);
            indexEF_store=indexEF_iterations; 
            Master_matrix_static(indexEF_store,1)=alpha1;
            Master_matrix_static(indexEF_store,2)=alpha2;
            Master_matrix_static(indexEF_store,3)=alpha3;
            Master_matrix_static(indexEF_store,4)=alpha4;
            Master_matrix_static(indexEF_store,5)=alpha5;
            Master_matrix_static(indexEF_store,6)=alpha6;
            Master_matrix_static(indexEF_store,7)=alpha7;
            %% GA Outputs
            Master_matrix_static(indexEF_store,8)=NaN;
            Master_matrix_static(indexEF_store,9)=NaN;
            Master_matrix_static(indexEF_store,10)=NaN;
            Master_matrix_static(indexEF_store,11)=NaN;
            %% Scaled Sector Totals (Static)
            Master_matrix_static(indexEF_store,12)=P1; 
            Master_matrix_static(indexEF_store,13)=P2;
            Master_matrix_static(indexEF_store,14)=P3;
            Master_matrix_static(indexEF_store,15)=P4;
            Master_matrix_static(indexEF_store,16)=P5; 
            Master_matrix_static(indexEF_store,17)=P6;
            Master_matrix_static(indexEF_store,18)=P7; 
            %% Sector Raw Totals 
            Master_matrix_static(indexEF_store,19)=mussel_NPV_tmp; % Mussel Raw 
            Master_matrix_static(indexEF_store,20)=Fish_NPV_tmp; % Finfish Raw 
            Master_matrix_static(indexEF_store,21)=kelp_NPV_tmp; % Kelp Raw 
            Master_matrix_static(indexEF_store,22)=Halibut_NPV_tmp; % Halibut Static Raw
            Master_matrix_static(indexEF_store,23)=Halibut_NPV_scaled_tmp; % Halibut Static Scaled 
            Master_matrix_static(indexEF_store,24)=Halibut_Yi_tmp; % Halibut Dynamic Raw 
            Master_matrix_static(indexEF_store,25)=halibut_dynamic_scaled_tmp; % Halibut Dynamic Scaled 
            Master_matrix_static(indexEF_store,26)=Viewshed_raw_value_fin_tmp; % Viewshed Finfish Raw 
            Master_matrix_static(indexEF_store,27)=Viewshed_raw_value_kelpmussel_tmp; % Viewshed Kelp/Mussel Raw 
            Master_matrix_static(indexEF_store,28)=enviro_impact_tmp; % Enviormental Impacts Raw 
            Master_matrix_static(indexEF_store,29)=disease_raw_sum;% Disease Raw  
            Master_matrix_static_itor=Master_matrix_static(indexEF_store,:);
%             cd('C:\Users\jsteve13\Desktop\Static MSP\Static Solutions')
%             eval(['save ', ['data_save_MSP_static_values_A1_', num2str(100*alpha1),'_A2_', num2str(100*alpha2),'_A3_', num2str(100*alpha3),'_A4_',num2str(100*alpha4),'_A5_', num2str(100*alpha5),'_A6_', num2str(100*alpha6),'_A7_', num2str(100*alpha7)],' Master_matrix_static_itor'])
%             cd('C:\Users\jsteve13\Desktop\Static MSP')
    end
    finish_out=struct('Store',Master_matrix_static);

