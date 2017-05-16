function [EXITFLAG_out, FVAL_out, generations_out, funccount_out,X_out]=...
    eval_parallel_funct3(alpha1,alpha2,alpha3,alpha4,alpha5,alpha6,alpha7)

%set directory of the data_save files
% STEM='C:\Users\jsteve13\Desktop\Evaluating Solutions\Evaluating Solutions'; %% Need to change depending on the directory 
% d=dir([STEM,'/','data_save_*']);
%set file name
% fname=[STEM,'/',d(1).,'A1_', num2str(100*alpha_1),'_A2_', num2str(100*alpha_2),'_A3_', num2str(100*alpha_3),'_A4_',num2str(100*alpha_4),'_A5_', num2str(100*alpha_5),'_A6_', num2str(100*alpha_6),'_A7_', num2str(100*alpha_7)];
fname=['data_save_MSP_static_values_A1_', num2str(100*alpha1),'_A2_', num2str(100*alpha2),'_A3_', num2str(100*alpha3),'_A4_',num2str(100*alpha4),'_A5_', num2str(100*alpha5),'_A6_', num2str(100*alpha6),'_A7_', num2str(100*alpha7)];
%load file
load(fname);

%output variables
EXITFLAG_out=EXITFLAG;
FVAL_out=FVAL;
if isnan(FVAL)==1
    generations_out=NaN;
    funccount_out=NaN;
else
generations_out=OUTPUT.generations;
funccount_out=OUTPUT.funccount;
end
X_out=X;

end