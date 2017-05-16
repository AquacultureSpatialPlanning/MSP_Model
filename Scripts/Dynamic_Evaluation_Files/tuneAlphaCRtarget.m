function [out] = tuneAlphaCRtarget(x)

global habitat_area_i numpatches max_age Bij age_mature Dii K age_legal Ei delta Nij Mii age_move Wij TargetBiomass beta_i...
    Rmax CRgoal Nij_initial theta price Fsum_msy_guess SSB_msy_target Theta_Density_prop Fsum alphaCR gamma distance_to_port_for_each_soft_depth_patch yield_scalar...
    Nij_initial_virgin CRtarget SQ_1fishable_0notfishable_for_each_soft_depth_patch qi omega nearestdist_to_outfall tol_simulation T x0_PPUE_initial discount_rate
 
alphaCR=x;
alphaCR(alphaCR<0)=0;
beta_i=alphaCR./(Rmax.*habitat_area_i); %determine beta

%run w/o harvest so can calculate theta (for Fsum opt in Rmax) and CR (for this opt)
Ei=zeros(numpatches,1); %no harvest
theta=0; %placeholder
BasicSpatialHarvestModel_COREv1_TS %run model to virgin K
PostHstockdensity_kvirgin_i=PostHstockdensity_i; %stock biomass density per patch
PostHstockdensity_kvirgin=(sum(PostHstockdensity_kvirgin_i.*habitat_area_i))/sum(habitat_area_i); %average across all habitat patches of stock biomass density per unit area habitat
theta=PostHstockdensity_kvirgin*Theta_Density_prop*price; %wrt LEGAL fish
CR;
out=(CR-CRtarget)^2;
end
