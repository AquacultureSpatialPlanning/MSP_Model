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
else
    Dynamic_values=csvread('Dynamic_values.csv');
end
