% %Crow White: This is a basic age-structured spatially explicit model that
% %can be used for a zillion different applications. I've commented it
% %thoroughly for clarity. For more explanation see Section 4.2 in the SI of 
% % White et al. 2012 PNAS. (Note: Different than Eq. S11 [Section 4.2.5] here
% % natural and fishing mortality are calculated using an exponential function). 
% %READ AND AGREE TO BEFORE CONTINUING: If you use this model for a
% %project/publication let me know. Also, provide me with the
% %opportunity to be a part of the project/publication as a co-PI or co-Author. Thanks! Crow

%Nomenclature of subscripts to parameters used below
% i = patch = rows
% j = age = columns

if T>499 %if running model to equilibrium
    Eiy=repmat(Ei,1,T); %then Ei is temporally constant, so set it so here
end

% %initial conditions
Nij=Nij_initial; %Initial number of fish. rows=patches; cols=size classes
Bij=Nij.*Wij; %biomass abundance in each patch and age class

%% simulation
%notes:
%1. The below functions represent one order of processes (Adults
%spawn-Adults harvested-Adults move). Other orders are of course possible
%2. It is important to think carefully (and keep track of) NUMBERS vs.
%DENSITIES. E.g., Density is important to recruitment whereas Number is
%important to fishery yield
% tic %start timer
%Set resolution of model soln representing equilibrium conditions:
% tol_simulation=1e-6; %proportion difference in biomass allowed. 1e-6 default
% tol_simulation=1e-5; %proportion difference in biomass allowed. testing...
%differences less than this (i.e., equilibrium has essentially been reached)
Kdiff=1; %Set dummy starting biomass difference
% max_simulation_length=1000; %Make a large number of years that won't be reached before tolerance is met

%create empty shells of the spatial and temporal data you wish to record - TIME SERIES DATA
Bijy=NaN(numpatches,max_age,T); %age-structured biomass
Yiy=NaN(numpatches,T); % yield
Payoffiy=NaN(numpatches,T); %payoff (e.g., profit)

%These calculations are constant across the simulation, so do now:
% Ei=Ei.*SQ_1fishable_0notfishable_for_each_soft_depth_patch; %set to zero effort in nonfishable patches
% Survivei_legal=exp(-(Ei+delta)); %survival rate of harvested fish in each patch ~ fishing + natural mortality. Same for all legal-to-fish ages
% Surviveij_legal=repmat(Survivei_legal,1,length(age_legal:max_age)); %survival rate per patch, per age class
% Survivei_nonlegal=exp(-(Ei.*0+delta)); %survival rate of non-harvested fish ~ natural mortality
% Surviveij_nonlegal=repmat(Survivei_nonlegal,1,length(1:(age_legal-1))); %survival rate per patch, per age class
% Surviveij=[Surviveij_nonlegal,Surviveij_legal]; %stitch survival rates together for all age classes    

y=0; %keeps track of years  - TIME SERIES DATA
for t=1:T %start of simulation loop - TIME SERIES DATA
    y=y+1; %increase by one each year - TIME SERIES DATA
    Bijy(:,:,y)=Bij; %biomass at beginning of year - TIME SERIES DATA
%     Kold=sum(Bij(:)); %store for comparison to end of this year (i.e., beginning of next year)
    
    Bi_mature=sum(Bij(:,age_mature:end),2); %biomass of mature fish in each patch
    Li=Bi_mature; %number of larvae produced in each patch ~ fish biomass
    Si=(Li'*Dii)'; %disperse the larvae according to transitionMeeting probabilities
    RNi=(Si.*alphaCR)./(1+Si.*beta_i); %number of recruits in each patch following Density Dependence
    %Note: the above DD function is Beverton-Holt (could be Ricker, etc.), and 
    %the beta_i includes area habitat area thus converts the Si in the denominator to *density* of settlers.
    
    Ei=Eiy(:,y); %pull spatial fishing effort for year t
    %calculate survival implications wrt natural and fishing mortality
    Survivei_legal=exp(-(Ei+delta)); %survival rate of harvested fish in each patch ~ fishing + natural mortality. Same for all legal-to-fish ages
    Surviveij_legal=repmat(Survivei_legal,1,length(age_legal:max_age)); %survival rate per patch, per age class
    Survivei_nonlegal=exp(-(Ei.*0+delta)); %survival rate of non-harvested fish ~ natural mortality
    Surviveij_nonlegal=repmat(Survivei_nonlegal,1,length(1:(age_legal-1))); %survival rate per patch, per age class
    Surviveij=[Surviveij_nonlegal,Surviveij_legal]; %stitch survival rates together for all age classes   
    
    %Biomass before harvest and natural mortality:
    Bi_legal=sum(Bij(:,age_legal:end),2); %LEGAL fish in each patch
%     Bi_ALL_preH=sum(Bij,2); %ALL size (set by age) fish in each patch
     
%     Yield - TIME SERIES DATA  
    Yi=(Bi_legal.*(1-Survivei_legal)).*(Ei./(Ei+delta)); %biomass yield in each patch - TIME SERIES DATA
    Yiy(:,y)=Yi; %biomass yield in each patch - TIME SERIES DATA

    %Escapement: biomass after harvest and natural mortality:
    Bij_escape=Bij(:,age_legal:end).*Surviveij(:,age_legal:end); %LEGAL fish in each patch
    % Bij_escape=Bij.*Surviveij; %ALL size (set by age) fish in each patch
    
    SSBij_postharvest=Bij(:,age_mature:end).*Surviveij(:,age_mature:end); %SSB right after harvest
    SSB_postharvest=sum(sum(SSBij_postharvest));
    
    %Surviving fish transition to next age class:
    Nij(:,2:end)=Nij(:,1:(end-1)).*Surviveij(:,1:(end-1)); %number of fish that survive to next age class (all fish at max age die)
    Nij(:,1)=RNi; %add recruit number to the first age class 
    % adult movement
    Nij_move=(Nij'*Mii)';
    Nij(:,age_move:max_age)=Nij_move(:,age_move:max_age); %update Nij with adults following their movement
    Bij=Nij.*Wij; %convert # to biomass
    K=sum(Bij(:)); %total biomass abundance
%     Kdiff=abs(K-Kold)/K; %see if equilibrium reached yet

    %Stock densities in each year
    PreHstockdensity_i=Bi_legal./habitat_area_i; %Biomass density of LEGAL fish before harvest and natural mortality
    PreHstockdensity_i(PreHstockdensity_i==Inf)=0; %set infinity patches (i.e., where habitat_area_i=0) to zero
%     PreHstockdensity_i=Bi_ALL_preH./habitat_area_i; %Biomass density of ALL fish before harvest (and natural mortality)
    PostHstockdensity_i=sum(Bij_escape,2)./habitat_area_i; %Biomass density of LEGAL fish after harvest AND natural mortality
    PostHstockdensity_i(PostHstockdensity_i==Inf)=0; %set infinity patches (i.e., where habitat_area_i=0) to zero
    %Economics
    Tcost_density_no_distance_effect=theta.*(log(PreHstockdensity_i)-log(PostHstockdensity_i)); %total cost per unit area (ie density)
    Tcost_density=Tcost_density_no_distance_effect.*(1+gamma.*distance_to_port_for_each_soft_depth_patch);
    Tcost=Tcost_density.*habitat_area_i; %convert total cost per unit area to total cost in each patch
    Tcost_fishery=Tcost.*(Ei./(Ei+delta)); %stock effect cost that fishery bears.  Q: IS THIS CORRECT? OR, SHOULD COST TO THE FISHERY ONLY BE A FCN OF EFFORT (NOT EFFORT TIMES CATCHABILITY)??? <<If change this then need to change fleet model function Ei_wrt_PPUE_qi
    Trevenue=Yi.*price;  
    %Calc and record profit
    Tprofit=Trevenue-Tcost_fishery; 
    Tprofit(PreHstockdensity_i==0)=0; %otherwise cost is infinity and profit =NaN
    Tprofit(PostHstockdensity_i==0)=0;%otherwise cost is infinity and profit =NaN
    Payoff=Tprofit;
    Payoffiy(:,y)=Payoff; %profit in each patch - TIME SERIES DATA
    % disp(['sum(Payoff)=',num2str(sum(Payoff))]);
    
end %end of simulation loop
  
%Calculate the compensation ratio (make sure Ei=0 for correct calculation)
% a=origin
% b=from origin to max number of settlers
%CR=(slope a)/(slope b)
% slope a (i.e, at S=0)
slope_a=alphaCR;
% slope b (i.e., at S=Svirgin=Smax)
slope_b=sum(RNi)/sum(Si);
%note: also works for slope_b=sum(RNi(n))/sum(Si(n));
CR=slope_a/slope_b;

%get rid of NaNs at end of records
% Yiy=Yiy(:,1:y); %- TIME SERIES DATA
% Bijy=Bijy(:,:,1:y);% - TIME SERIES DATA
% Payoffiy=Payoffiy(:,y);  %- TIME SERIES DATA
% disp(['Simulation took ',num2str(toc),' seconds'])

% plot_TS_results