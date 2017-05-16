% load('Tuner_save_Jun_22_2015_15_28','Aqua_dev_indices','mussel_NPVi','Fish_NPVi','kelp_NPVi','Halibut_10yrNPVi','Viewshed_max_impact','enviro_impact','V2','disease_metrici')
index_manual=1:25;
% % Number of Jobs 
%     number_jobs=25;
% % Number of Tasks 
%     number_tasks=25;
% % Number of 125
%     number_iterations=125;
% % Create Compiler 
% Static_compiler_matrix=NaN(25,25,125);
% itor_count1=0;
% itor_counter=0;
% while itor_count1<size(Static_compiler_matrix,1)
%     itor_count1=itor_count1+1;
%     itor_count2=0;
%     while itor_count2<size(Static_compiler_matrix,2)
%         itor_count2=itor_count2+1;
%         itor_count3=0;
%         while itor_count3<size(Static_compiler_matrix,3)
%             itor_counter=itor_counter+1;
%             itor_count3=itor_count3+1;
%             Static_compiler_matrix(itor_count1,itor_count2,itor_count3)=itor_counter;
%         end
%     end
% end % Clear in order to increase computational efficiancy 
% load('Static_MSP_solutions')
% disp('Starting Jobs......')
c=parcluster('local');
for index_jobs=index_manual%:size(Static_compiler_matrix,1)
    jfinish=createJob(c); 
    for index_tasks=1:size(Static_compiler_matrix,2)
        index_iterations=1:size(Static_compiler_matrix,3);
        indexEF_tasks_compile=Static_compiler_matrix(index_jobs,index_tasks,index_iterations);
        indexEF_tasks=indexEF_tasks_compile(:); X_in=EVwschoice(:,indexEF_tasks);
%         vector_indexes=[index_jobs;index_tasks;index_iterations];
        field1='Task_solutions';value1=X_in;field2='Alphas';value2=aMatrix;field3='Job_index';value3=index_jobs;
        field4='Task_index';value4=index_tasks;field5='Iteration_index';value5=index_iterations;field6='Index_mat';value6=indexEF_tasks;
        tasks_in=struct(field1,value1,field2,value2,field3,value3,field4,value4,field5,value5,field6,value6);
%         Finishing_function_static(tasks_in)
        t=createTask(jfinish,@Finishing_function_static, 1, {tasks_in}, 'CaptureDiary',true); 
%         finishTime=toc(timerVal);
%         disp(num2str(finishTime))
%         disp(['Iteration set ',num2str(indindex_jobs),' submitted'])
%         finishTime=toc(timerVal);
%         disp(['took ',num2str(finishTime),' seconds'])
    end
    % finishTime=toc(timerVal)
    % finishTime=toc(timerVal);
    % disp(['took ',num2str(finishTime),' seconds'])
    jfinish.submit;
%     disp(num2str(finishTime))
    disp(['Iteration set ',num2str(index_jobs),' submitted'])
end
% disp('Job Finished...........')
