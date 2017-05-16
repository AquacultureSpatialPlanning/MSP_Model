    %% Finishing Function
    function finish_out=Finishing_function_static(tasks_in)
    input_data_dir = tasks_in.input_data_dir;
    % load(strcat(input_data_dir,'tuned_params'))
    load(strcat(input_data_dir,'Aqua_Dev_Indices'))
    load('tuned_parameters_short')
    warning('off','all')
    gammatmp=gamma_var;
    % cd(strcat(script_dir,'Dynamic_Evaluation_Files/'))
    % addpath(input_data_dir)
    % Load all data for the task specific structure
    X_iter=tasks_in.Task_solutions;
    indexEF_tasks=tasks_in.Index_mat;
    Halibut_Dynamic=NaN(length(indexEF_tasks),1);
    for indexEF_iterations=1:length(tasks_in.Iteration_index)
            X=X_iter(:,indexEF_iterations);
            indexEF=indexEF_tasks(indexEF_iterations);
            [Halibut_Yi_tmp]=Sector_calcs_wrt_X(X,input_data_dir,target_fid_fulldomain,Aqua_Dev_Indices,SQ_1fishable_0notfishable_for_each_soft_depth_patch_ORIGINAL,Nij_initial_msy_TS,x0_PPUE_msy_TS,Fsum_msy,target_fid_hab_depth,Wij,numpatches,max_age,T,delta,age_legal,age_mature,Dii,alphaCR,beta_i,habitat_area_i,theta,price,gammatmp,distance_to_port_for_each_soft_depth_patch,Mii,age_move,discount_rate_iy);
            indexEF_store=indexEF_iterations;
            Halibut_Dynamic(indexEF_store,1)=Halibut_Yi_tmp; % Halibut Dynamic Raw
            % Master_matrix_static_itor=Master_matrix_static(indexEF_store,:);
    end
    finish_out=struct('Store',Halibut_Dynamic);
