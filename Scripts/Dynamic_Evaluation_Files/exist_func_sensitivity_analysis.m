function [exist_var,exist_vec]=exist_func_sensitivity_analysis(exist_vec,indexEF,alpha1,alpha2,alpha3,alpha4,alpha5,alpha6,alpha7,solution_directory,script_directory,solution_name)

%set directory of the data_save files
% STEM='C:\Users\jsteve13\Desktop\Evaluating Solutions\Evaluating Solutions'; %% Need to change depending on the directory 
% d=dir([STEM,'/','data_save_*']);
%set file name
% fname=[STEM,'/',d(1).,'A1_', num2str(100*alpha_1),'_A2_', num2str(100*alpha_2),'_A3_', num2str(100*alpha_3),'_A4_',num2str(100*alpha_4),'_A5_', num2str(100*alpha_5),'_A6_', num2str(100*alpha_6),'_A7_', num2str(100*alpha_7)];
%fname=['data_save_A1_', num2str(100*alpha1),'_A2_', num2str(100*alpha2),'_A3_', num2str(100*alpha3),'_A4_',num2str(100*alpha4),'_A5_', num2str(100*alpha5),'_A6_', num2str(100*alpha6),'_A7_', num2str(100*alpha7)];
cd(solution_directory)
%load file
if (exist([[solution_name, num2str(100*alpha1),'_A2_', num2str(100*alpha2),'_A3_', num2str(100*alpha3),'_A4_',num2str(100*alpha4),'_A5_', num2str(100*alpha5),'_A6_', num2str(100*alpha6),'_A7_', num2str(100*alpha7)],'.mat'],'file'))==0 
    exist_var=0;
else 
    exist_var=1;
end
cd(script_directory)
exist_vec(indexEF)=exist_var;