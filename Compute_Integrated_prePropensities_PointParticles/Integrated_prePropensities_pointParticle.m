function [Integrated_CDFs,distances,numVariances,var_list] = Integrated_prePropensities_pointParticle(...
    contact_sphere_radiusSqd,T_HARDLIMIT_DIFFUSE,maxRxDist_T_HARD,Da,Db,...
    earliestPropensityTime)

% Point Particle Representation
% A prePropensity is just the propensity/rateConstant

% % COMPUTE INPUTS
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

% % FIXED PARAMETERS - old
% n_sigma = 5; % used to calc max diffusion distance corresponding to max diffusion time
% Dfast = 1; 
% contact_sphere_radiusSqd = .8;
% T_HARDLIMIT_DIFFUSE = 5; % units [s]
% Da = 1; Db = 1; % units [um^2 / s]
% maxRxDist_T_HARD = n_sigma*sqrt(6*Dfast*T_HARDLIMIT_DIFFUSE); % units [um]
% % FIXED PARAMETERS

numDistances = 2000;%800;
numVariances = 5000;%1600;

d_min = sqrt(contact_sphere_radiusSqd);
d_max = 2*maxRxDist_T_HARD; % each particle can diffuse a max distance <maxRxDist_T_HARD>. Thus max separation dist is twice that.

d_list = linspace(d_min,d_max,numDistances);
distances = d_list;

d_list = repmat(d_list,numVariances,1);
d_list = d_list(:);

% Used for numerical integration of pre-propensities
times_list = linspace(earliestPropensityTime,T_HARDLIMIT_DIFFUSE,numVariances)';   


var_list = (6*Da*times_list) + (6*Db*times_list);

% Outputs: both are matrices. Each row corresponds to a pre-propensity  
% (/integrated pre-propensity) curve at a given distance.
[CDFs,Integrated_CDFs] = compute_Integrated_prePropensity_curves(...
     contact_sphere_radiusSqd,distances,var_list,times_list);
 save(['SAVED_propensities__and_prePropensities_pointParticle_',date,'.mat'])
 
end

