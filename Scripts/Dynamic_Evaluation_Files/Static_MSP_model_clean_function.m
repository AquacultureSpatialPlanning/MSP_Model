% close all
% clear all
function Static_MSP_model_clean_function
%this code show in concept how to conduct a 7D tradeoff analysis.
%the variables are modeled after the SCB aquaculture tradeoff study, but
%all the values are fake numbers
% prompt='Epsilon = ';
% str=input(prompt,'s');
epsilon=.2;
tic
load('Tuner_save_Jun_22_2015_15_28','Aqua_dev_indices','mussel_NPVi','Fish_NPVi','kelp_NPVi','Halibut_10yrNPVi','Viewshed_max_impact','enviro_impact','V2','disease_metrici')
N=1061; %number of patches that can be developed for aqua - Can change this
nS=7; %number of sectors
%Raw values of each sector in each patch
% Emerging sectors
M=mussel_NPVi(Aqua_dev_indices); % 90% of patches can be developed for M
F=Fish_NPVi(Aqua_dev_indices); % 30% of the patches can be developed for F, some are outside M
K=kelp_NPVi(Aqua_dev_indices); % 10% of the patches can be developed for K. F and K do not overlap
% Existing sectors
H = Halibut_10yrNPVi(Aqua_dev_indices); %Value that would be lost if M, F or K were developed in the patch
V = Viewshed_max_impact(Aqua_dev_indices); %Impact in patch if M, F or K were developed there
B = enviro_impact(Aqua_dev_indices); %Impact that would happen if F was developed there
D=zeros(length(B),1);D(V2(Aqua_dev_indices))=disease_metrici; %Impact that would happen if F was developed there
%Matrix of the potential raw value of each sector in each patch
EVr=[M, F, K, H, V, B, D];
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
Vs=V./sum(V);
Bs=B; Bs(Fs==0)=0; Bs=Bs./nansum(Bs); %first remove values that can't be affected
Ds=D; Ds(Fs==0)=0; Ds=Ds./nansum(Ds);%first remove values that can't be affected

%Set sector weighting scenarios
disp('Set sector weighting scenarios')
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
    aKi+aKi;
    aHi=0;
    for aH=a_range
        aHi=aHi+1;
        aVi=0;
        for aV=a_range
            aVi+aVi;
            aBi=0;
            for aB=a_range
                aBi+aBi;
                aDi=0;
                for aD=a_range
                    aDi+1;
                    iMatrix=iMatrix+1;
                    aMatrix(iMatrix,:)=[aM,aF,aK,aH,aV,aB,aD];
                end
            end
        end
    end
end
end
end

%Solve the sector weighting scenarios
disp('Solve the sector weighting scenarios')
shell=NaN(N,length(aMatrix(:,1)));
Policy=shell;
EVwschoice=shell;
EVr_wrt_N_alpha=NaN(N,nS,length(aMatrix(:,1)));
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
EVsw2 = Ms.*aM + Bs.*aB + Ds.*aD; %Develop M
EVsw3 = Fs.*aF; %Develop F
EVsw4 = Ks.*aK + Bs.*aB + Ds.*aD; %Develop K

%Identify which policy is best for each patch
EVswoptions=[EVsw1 EVsw2 EVsw3 EVsw4];
[Y, I] = max(EVswoptions,[],2); % I indicates which policy to implement in each patch, given the weights
Policy(:,ai)=I; %rows=patch policy number; col=weighting scenario
EVwschoice(:,ai)=Y; %rows=patch policy weighted scaled value; col=weighting scenario
end
disp(['FINISHED: Took ',num2str(toc/60/60),' hours'])
% Convert policy integers to ones which make more intuitive sense (i.e.
% convert (1 to 0 for no development, etc)
Static_plans=Policy;
Static_plans(Policy==1)=0;Static_plans(Policy==2)=1;Static_plans(Policy==3)=2;Static_plans(Policy==4)=3;
Static_values=NaN(size(Static_plans,2),7);
% Find the values of each sector for each developed policy
for index=1:size(EVwschoice,2);
    % Mussel 
        Static_values(index,1)=sum(Ms(Static_plans(:,index)==1));
    % Finfish
        Static_values(index,2)=sum(Fs(Static_plans(:,index)==2));
    % Kelp
        Static_values(index,3)=sum(Ks(Static_plans(:,index)==3));
    % Halibut 
        Static_values(index,4)=sum(Hs(Static_plans(:,index)==0));
    % Viewshed
        Static_values(index,5)=sum(Vs(Static_plans(:,index)==0));
    % Benthic 
        Static_values(index,6)=sum(Bs(Static_plans(:,index)~=2));
    % Disease 
        Static_values(index,7)=sum(Ds(Static_plans(:,index)~=2));
end
% save('Static_MSP_solutions','EVwschoice','Policy','Static_plans','aMatrix','Static_values')
%% Export Data to R 
    csvwrite('Static_MSP_plans.csv',Static_plans)
    csvwrite('aMatrix.csv',aMatrix)
    csvwrite('Policy.csv',Policy)
    csvwrite('Static_values.csv',Static_values)







