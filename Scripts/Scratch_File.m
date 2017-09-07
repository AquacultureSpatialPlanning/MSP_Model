% Figure 3
load ~/MSP_Model/Input/Data/Lester_et_al_MSPsolutions_Evals_v3/Policy_i_a.mat
Policy_i_a_A = Policy_i_a - 1; Policy_i_a_M = Policy_i_a - 1; Policy_i_a_F = Policy_i_a - 1; Policy_i_a_K = Policy_i_a - 1;
% Aquaculture
Policy_i_a_A(Policy_i_a_A ~= 0) = 1;
Freq_A = sum(Policy_i_a_A,2) / size(Policy_i_a_A,2);
% Mussel
Policy_i_a_M(Policy_i_a_M ~= 1) = 0;
Freq_M = sum(Policy_i_a_M,2) / size(Policy_i_a_M,2);
% Finfish
Policy_i_a_F(Policy_i_a_F ~= 2) = 0;
Policy_i_a_F(Policy_i_a_F == 2) = 1;
Freq_F = sum(Policy_i_a_F,2) / size(Policy_i_a_F,2);
% Kelp
Policy_i_a_K(Policy_i_a_K ~= 3) = 0;
Policy_i_a_K(Policy_i_a_K == 3) = 1;
Freq_K = sum(Policy_i_a_K,2) / size(Policy_i_a_K,2);

save('Fig3_data.mat','Freq_A', 'Freq_M', 'Freq_F', 'Freq_K')

% Figure 4
clear all
load ~/MSP_Model/Input/Data/Lester_et_al_MSPsolutions_Evals_v3/Policy_i_a_ALL_wrt_DM_f01.mat
Policy_i_CS = Policy_i_a_ALL_wrt_DM_f01 - 1;
Policy_i_CS(Policy_i_CS ~= 0) = 1;
Freq_CS = sum(Policy_i_CS,2) / size(Policy_i_CS,2);

save('Fig4_data.mat','Freq_CS')

% Figure 5
