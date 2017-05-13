function [out]=habitat_area_i_wrt_habitat_prop_reduction_wrt_depth(psi,depth_i,habitat_area_i_raw)

%use psi in a negative exponential function modeled after Haugen 1990
% Parameters based on fitting a log function to the empirical data. Params: psi=-0.233; intercept=1.0384; R2=0.9658
habitat_prop_reduction_wrt_depth=psi.*log(-depth_i+1)+1.0384;
%get rid of negative values
habitat_prop_reduction_wrt_depth(habitat_prop_reduction_wrt_depth<=0)=0; %Set no zero.
habitat_area_i=habitat_area_i_raw.*habitat_prop_reduction_wrt_depth; %modify habitat area wrt the reduction function
out=habitat_area_i;