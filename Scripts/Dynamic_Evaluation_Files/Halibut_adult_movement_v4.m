% global numpatches habitat_area_i

disp('Calculate adult movement')

    %Normalized Gaussian function of 2D probability distribution wrt to
    %Diffusivity (D), elapsed time (t), and distance from centroid (r). 
    %See below for a solution to the 2D isotropic diffusion equation:
    %see: http://www.google.com/url?sa=t&rct=j&q=&esrc=s&source=web&cd=1&ved=0CDUQFjAA&url=http%3A%2F%2Fwww.rpgroup.caltech.edu%2Fcourses%2Faph162%2F2006%2FProtocols%2Fdiffusion.pdf&ei=2MJIUZmbNKSAiwLyp4DQAg&usg=AFQjCNF2G6VnVLeMDQ_s9EsROoSMfMsAJg&bvm=bv.44011176,d.cGE
    %Above document is saved as diffusion-1.pdf in H:\Projects\EEZ and MPAs\Global\Global_Grid_Logisticv1
    %Note from Tristan:
    %    The formula in the original document is identical to the Green's function in the newer document ... 
    %    it's just getting rid of the integration variable x' ... i.e. set x'=0, and then |x|^2 = x^2 + y^2.  
    %    Also, the other document you sent has the same form, except they use alpha^2=D, and it is for 1D not 2D, 
    %    which gives the square root in the denominator (see eqn 60 in the newer tutorial).  Also, the t' in 
    %    that other presentation is a dummy variable, which generally just confuses things, 
    %    but it can be set to zero, i.e. t'=0, and t>=0.
    % The 2D isotropic diffusion equation is known as "Green's function"
    %References:
    % http://mathworld.wolfram.com/GreensFunction.html
    % Arfken, G. "Nonhomogeneous Equation--Green's Function," "Green's Functions--One Dimension," and "Green's Functions--Two and Three Dimensions." §8.7 and §16.5-16.6 in Mathematical Methods for Physicists, 3rd ed. Orlando, FL: Academic Press, pp. 480-491 and 897-924, 1985.
    % Garabedian, P. R. Partial Differential Equations. New York: Wiley, 1964.
    % Marichev, O. I. "Funktionen vom hypergeometrischen Typ und einige Anwendungen auf Integral- und Differentialgleichungen." Ph.D. dissertation. Jena, Germany: Friedrich-Schiller-Universität, p. 266, 1990. 
    % http://en.wikipedia.org/wiki/Green%27s_function
%% CA halibut adult movement estimate wrt empirical estimate of daily migration rate
    
%Halibut demographics
D=0.21; %km/day (Domeier and Chun 1995 CalHal_movement)
t=365; %days per year
k=2; %k is a scaling factor representing the sensitivity of the calculated emigration rate to changes in environmental
% suitability (as measured by D) (Walters et al. 1999). Small values of k (e.g. 0.1) result in high
% sensitivity to change in D while large values (e.g. 10) render adult dispersal rate insensitive to D.
% We used an intermediate value of k (= 2), but we also compared model outputs from alternative k
% values. See Cheung2008_ModelingFishDistributionReport

%data shells
shell=NaN(numpatches,numpatches);
r_patch_pairwisekmdistance=shell;
G_patch_pairwiseGreensfcn=shell;
Hr_patch_pairwise_hab_suitability=shell;
mij=shell;

%For every focal cell...
for a=1:numpatches
%     disp(['Patches left: ',num2str(numpatches-a)])
    %calculate distance from it to all patches in the domain...
    r=deg2km(distance(lat_lon_habitat_patches(a,:),lat_lon_habitat_patches)); r_patch_pairwisekmdistance(a,:)=r; %pairwise distance between cells in km
    G=exp((-r.^2)./(4*D*t))./(4*pi*D*t); %green's function
    G_rowsums2one=G./sum(G); %rowsum standardized so that all the fractions of movement 
%     from the focal cell to the other cells sum to ONE. Doing this before
%     introducing habitat suitability effects so that the Green's fcns
%     values are not tiny and overwhelmed by the D values
    G_patch_pairwiseGreensfcn(a,:)=G_rowsums2one; 
    Hr=habitat_area_i(a)./habitat_area_i; 
    Hr_patch_pairwise_hab_suitability(a,:)=Hr; %relative habitat area
    mij_raw=(G_rowsums2one.*k)./(k+Hr); %raw migration rates from focal patch to all patches
    mij(a,:)=mij_raw./sum(mij_raw); %rowsum standardized so that all the fractions of movement 
%     from the focal cell to the other cells sum to ONE
end

%plot and save data
figure; 
surface(mij); 
title('Adult movement'); 
xlabel('Destination'); 
ylabel('Source'); 
set(gcf,'color',[1 1 1]); 
colorbar; 
shading flat; 
axis tight
axis square
savefig('m_i_j')

%plot a particular focal patch results
a=2000;
figure
plot(r_patch_pairwisekmdistance(a,:),mij(a,:),'b.')
xlabel('Distance from focal patch [km]')
ylabel('m_i_j')
set(gcf,'color','white'); 
box off
title('i = source; j = destination')
savefig('Adult_movement_wrt_pairwisepatchdistance')

figure
semilogx(Hr_patch_pairwise_hab_suitability(a,:),mij(a,:),'b.')
xlabel('Relative habitat suitability: H_i / H_j')
ylabel('m_i_j')
set(gcf,'color','white'); 
box off
title('i = source; j = destination')
savefig('Adult_movement_wrt_pairwisepatchhabitatsuitability')

figure
plot((k./(k+Hr_patch_pairwise_hab_suitability(a,:))),mij(a,:),'b.')
xlabel('k / (k + H_i / H_j)')
ylabel('m_i_j')
set(gcf,'color','white'); 
box off
% xlim([0 10])
title('i = source; j = destination')
savefig('Adult_movement_wrt_pairwisepatchhabitatsuitabilityfcn')

% data_save_adult_movement=mij;
% save data_save_adult_movement
if halibut_max_depth==-183
    save('mij_183depthlimit','mij')
elseif halibut_max_depth==-90
    save('mij_90depthlimit','mij')
else
    save('mij','mij')
end

%%End of halibut code
%for more code on estimating movement using PDFs see the EEZ high seas
%study wrt grid cell locations!
