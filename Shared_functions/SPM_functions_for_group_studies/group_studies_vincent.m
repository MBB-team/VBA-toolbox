%%%%%%%%%%%%%%% COMPARISON - RFX RANDOM EFFECTS

% SPM Files required
%
% - spm_BMS_gibbs
% FORMAT [exp_r,xp,r_samp,g_post] = spm_BMS_gibbs (lme, alpha0, Nsamp)
%
%  INPUT:
%  lme      - array of log model evidences 
%               rows: subjects
%               columns: models (1..Nk)
%  alpha0   - [1 x Nk] vector of prior model counts
%  Nsamp    - number of samples (default: 1e6)
%  
%  OUTPUT:
%  exp_r   - [1 x  Nk] expectation of the posterior p(r|y)
%  xp      - exceedance probabilities
%  r_samp  - [Nsamp x Nk] matrix of samples from posterior
%  g_post  - [Ni x Nk] matrix of posterior probabilities with 
%            g_post(i,k) being post prob that subj i used model 

% - spm_gamrnd  ( A sampling function ) 
% - spm_multrnd  ( A sampling function ) 

%%%%%%%%%%%%%%% COMPARISON - FFX FIXED EFFECTS

% - spm_api_bmc (main function called)
% function out = spm_run_bms_dcm (varargin)
% API to compare DCMs on the basis of their log-evidences. Four methods
% are available to identify the best among alternative models:
%
%  (1) single subject BMS using Bayes factors
%     (see Penny et al, NeuroImage, 2004)
%  (2) fixed effects group BMS using group Bayes factors
%     (see Stephan et al, NeuroImage, 2007)
%  (3) random effects group BMS using exceedance probabilities
%     (see Stephan et al, NeuroImage, 2009)
%  (4) comparing model families
%     (see Penny et al, PLOS-CB, submitted)
%
% Note: All functions use the negative free energy (F) as an approximation
% to the log model evidence.


% - smp_api_bmc (function called when family comparison - see literature)
