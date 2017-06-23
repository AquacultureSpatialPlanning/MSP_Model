
  %% Evaluate dynamic response to each MSP aquaculture solution
  p = mfilename('fullpath');
  part = fileparts(p);
  parts = strsplit(part, '/Scripts/Halibut');
  DirPart = parts{end-1};
  tic
  disp('Evaluate MSP solutions')
  addpath(genpath(DirPart))
  load(strcat(DirPart,'/Input/Data/Raw_Impacts_FID.mat'))
  load(strcat(DirPart,'/Input/Data/tuned_params.mat'))
  %1061x1 vector of FIDs for each aqua site
  FID_1061aqua=double(Raw_Impacts.FID); %FIDs of the 1061 aqua cells

  %FIDs of the 4552 with halibut
  target_fid_fulldomain_wrt_filter_habitat_rows=target_fid_fulldomain(filter_habitat_rows==1);

  %Determine FIDs in common between halibut model and aqua models
  [C2,IA2,IB2] = intersect(target_fid_fulldomain_wrt_filter_habitat_rows,FID_1061aqua,'stable');
  % IA2 = 1031 rows of 4552x1 target_fid_fulldomain_wrt_filter_habitat_rows with FIDs matching 1061x1 FID_1061aqua
  % IB2 = 1031 rows of 1061x1 FID_1061aqua with FIDs matching 4552x1 target_fid_fulldomain_wrt_filter_habitat_rows

  %Load MSP solutions (279936)
  %load Policy_i_a.mat %1=ND; Greater than 1 = aqua devel
  Policy_i_a_tmp = load(strcat(DirPart,'/Input/Data/C_C_obj_i.mat'));
  Policy_i_a = Policy_i_a_tmp.C_C_obj_i + 1;%struct2cell(
  % Policy_i_a = NaN(1061,1062);
  % for itor = 1:size(Policy_i_a_tmp,1)
  %    Policy_i_a(:,itor) = Policy_i_a_tmp{itor,1} + 1;
  % end
  %Find unique policies (48268 out of 279936)
  Policy_i_a_trans=Policy_i_a';
  [Policy_i_a_trans_C,Policy_i_a_trans_IA,Policy_i_a_trans_IC] = unique(Policy_i_a_trans,'rows','stable');
  % Policy_i_a_trans_C          48268x1061  Indicates the unique values in A
  % Policy_i_a_trans_IA         48268x1    Indicates the first row in A where each C occurs
  % Policy_i_a_trans_IC        279936x1  Indicates what rows in C to fill each row in A
  % C = A(IA,:) and A = C(IC,:)

  %Choose below params
  T=10; %number of years to run (T=10)
  discount_rate=0.05; %discount rate (d=0.05)
  discount_factor=1./((1+discount_rate).^[1:T]);
  discount_rate_iy=repmat(discount_factor,numpatches,1);

  Y_NPV_wrt_MSP=NaN(length(Policy_i_a(1,:)),1); %shell
  Y_NPV_wrt_MSP_unique=NaN(length(Policy_i_a_trans_C(:,1)),1); %shell
  ipt_counter=round(linspace(1,length(Policy_i_a_trans_C(:,1)),10));
  warning ('off','all');
  for pt=1:length(Policy_i_a_trans_C(:,1))
    disp(pt)
%       if isempty(intersect(pt,ipt_counter))==0
%           disp(['pt = ',num2str(pt)]) %counter
%           save DynHalEval_wrt_MSP Y_NPV_wrt_MSP Y_NPV_wrt_MSP_unique pt %save data generated so far
%       end

      %Index the unique policy and transpose it
      Policy_i_a_unique=(Policy_i_a_trans_C(pt,:))';

      %Set unfishable patches wrt aqua policy
      SQ_1fishable_0notfishable_for_each_soft_depth_patch=SQ_1fishable_0notfishable_for_each_soft_depth_patch_ORIGINAL;
      SQ_1fishable_0notfishable_for_each_soft_depth_patch(IA2(Policy_i_a_unique(IB2)>1))=0; %Note:Policy_i_a=1 means No Development

      Nij_initial=Nij_initial_msy_TS; %initial stock status (patch and age class specific)
      %Run with status quo MSY management and yes Aquaculture
      Fsum=Fsum_msy; %TAE (fleet model)
      x0_PPUE_initial=x0_PPUE_msy_TS; %initial guess for fleet model
      BasicSpatialHarvestModel_COREv1_fleet_TS  %Run simulations
      tmp1=Yiy.*discount_rate_iy; %multiply year and patch specific values by discounted value = NPV_i
      Yi_NPV_MSP=sum(tmp1,2); %sum years together to calc NPV for each patch
      Y_NPV_wrt_MSP_unique(pt)=sum(Yi_NPV_MSP); %sum across years and save unique policy result
  end
  save(strcat(DirPart,'/Output/Data/CC_Dynamic_Halibut.mat'),'Y_NPV_wrt_MSP_unique')
exit;
