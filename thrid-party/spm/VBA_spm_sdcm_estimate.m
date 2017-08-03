function [DCM] = VBA_spm_sdcm_estimate(P,stochastic,augment,DisplayWin,confounds)
% Estimate parameters of a stochastic DCM for fMRI data
% FORMAT [DCM] = spm_sdcm_estimate(DCM,stochastic,augment)
%
% DCM           - the DCM or its filename
% stochastic    - flag for s/d DCM
% augment       - structure for input basis function set (see dcm2vba.m)
% DisplayWin    - flag for inversion report window
% confounds     - onfounds structure
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Jean Daunizeau

try; stochastic; catch, stochastic=0; end
try; augment; catch, augment=0; end
try; DisplayWin; catch, DisplayWin=0; end
try; confounds; catch, confounds='find'; end

% load DCM structure
if ~nargin
    P = VBA_spm_select(1,'^DCM.*\.mat$','select DCM_???.mat');
    if isempty(P)
        return
    end
end
if isstruct(P)
    DCM = P;
    P = ['DCM-' date];
else
    load(P)
end

% export DCM structure to VBA inversion inputs
[y,u,f_fname,g_fname,dim,options] = dcm2vba(DCM,stochastic,augment,confounds);

disp(' ')
disp('-----------------------------')
disp('Inverting DCM...')

% invert stochastic DCM using discrete-time VB routine:
options.DisplayWin = DisplayWin;
[posterior,out] = VBA_NLStateSpaceModel(y,u,f_fname,g_fname,dim,options);

% extract relevant info and fills in DCM structure
[DCM] = vba2dcm(posterior,out,DCM);
% spm_dcm_explore(DCM)

disp('Inverting DCM... OK.')
disp('-----------------------------')
disp(' ')

return








