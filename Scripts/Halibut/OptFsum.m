function [out] = OptFsum(x)

global alphaCR habitat_area_i numpatches max_age Bij age_mature Dii K age_legal Ei delta Nij Mii age_move Wij TargetBiomass beta_i...
    Nij_initial theta price Theta_Density_prop K_legal_virgin_abundance PostHstockdensity_i PostHstockdensity_kvirgin...
    Kvirgin_abundance Fsum habitat_area_i Bi_legal Yi yield_scalar gamma distance_to_port_for_each_soft_depth_patch Nij_initial_virgin...
    SQ_1fishable_0notfishable_for_each_soft_depth_patch qi omega nearestdist_to_outfall tol_simulation T x0_PPUE_initial discount_rate
 
% global Fsum Nij_initial


    Fsum=x;
    Fsum(Fsum<0)=0;
%     disp(['Fsum=',num2str(Fsum),]) 
%Correct way
    BasicSpatialHarvestModel_COREv1_fleet_TS
%shortcut
% Ei=(Fsum/numpatches).*ones(numpatches,1);
% BasicSpatialHarvestModel_COREv1

    
%     out=-sum(Payoff);
% disp(['Fsum=',num2str(Fsum),'; K=',num2str(K),'; sum(Yi)=',num2str(sum(Yi)),])
        out=-sum(Yi);
end

% Fsum_range=linspace(800,1200,10);
% sumYi=NaN(size(Fsum_range));
% for index=1:length(Fsum_range)
%     Fsum=Fsum_range(index)
%     BasicSpatialHarvestModel_COREv1_fleet
%     sumYi(index)=sum(Yi);
% end
% figure; plot(Fsum_range,sumYi)


