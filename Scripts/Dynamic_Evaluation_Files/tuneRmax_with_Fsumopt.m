function [out] = tuneRmax_with_Fsumopt(x)

global alphaCR habitat_area_i numpatches max_age Bij age_mature Dii K age_legal Ei delta Nij Mii age_move Wij TargetBiomass beta_i...
    Nij_initial theta price gamma distance_to_port_for_each_soft_depth_patch yield_scalar Fsum_msy_guess Fsum Nij_initial_virgin Theta_Density_prop...
    Bi_legal actual_comm_plus_rec_kglandings Fsum_UB SQ_1fishable_0notfishable_for_each_soft_depth_patch qi Nij_Fsumchoice omega nearestdist_to_outfall...
    tol_simulation T x0_PPUE_initial discount_rate
 
Rmax=x;
beta_i=alphaCR./(Rmax.*habitat_area_i); %set beta

%first run w/o harvest so can set theta
Ei=zeros(numpatches,1); %set no harvest to can determine PostHstockdensity_kvirgin for calculating theta
Nij_initial=Nij_initial_virgin; %set initial conditions 
theta=0; %placeholder
BasicSpatialHarvestModel_COREv1_TS %run model to virgin K
Nij_initial_virgin=Nij; %reset Nij_initial_virgin (note: did this in psi too)
PostHstockdensity_kvirgin_i=PostHstockdensity_i; %stock biomass density per patch
PostHstockdensity_kvirgin=(sum(PostHstockdensity_kvirgin_i.*habitat_area_i))/sum(habitat_area_i); %average across all habitat patches of stock biomass density per unit area habitat
theta=PostHstockdensity_kvirgin*Theta_Density_prop*price; %set theta 
CR_wrt_Rmax=CR;       
        
%Find MSY with fleet model
Nij_initial=Nij_initial_virgin.*24; %initial guess
TolX_value=1e-4*numpatches; %Set tolerance for Fsum search. 1e-4 is the default TolX. 
TolFun_value=TolX_value;
% We have so many patches and thus Fsum is so large (and fracional changes in Fsum don't matter much) so let's make Fsum larger
options=optimset('Display','off','TolX',TolX_value,'TolFun',TolFun_value);
Fsum_UB=sum(habitat_area_i)*4.5e-6;
[Fsum_opt negprofit flag]=fmincon(@OptFsum,Fsum_msy_guess,[],[],[],[],0,Fsum_UB,[],options);
Fsum_msy_guess=Fsum_opt;       
        
%Now run with Fsum=Fsum_msy so can calculate yield
Fsum=Fsum_opt;
BasicSpatialHarvestModel_COREv1_fleet_TS
Nij_Fsumchoice=Nij; % reset Nij_Fsumchoice (also reset in psi tuner)
        
%Compare with actual data
actual_comm_plus_rec_kglandings=177779; %Actual CDFW data on kg landings of CA halibut in SoCal (see Maunder et al. 2011 CA halibut stock assessment figures.xlsx sheet "Fig S2 + Fig S4 Other")
disp(['Given Rmax try = ',num2str(Rmax),' and alphaCR=',num2str(alphaCR),'; theta=',num2str(theta),'; CR_wrt_Rmax=',num2str(CR_wrt_Rmax),'; at Fsum_opt=',num2str(Fsum),' actual_comm_plus_rec_kglandings-sum(Yi) = ',num2str(actual_comm_plus_rec_kglandings-sum(Yi)),])
out=(actual_comm_plus_rec_kglandings-sum(Yi))^2; %for fmin algorthm
end
