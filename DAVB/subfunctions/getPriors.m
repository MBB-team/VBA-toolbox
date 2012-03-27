function [priors] = getPriors(nreg,n_t,options,reduced_f,stochastic)
% builds priors for DCM inversion
% function [priors] = getPriors(nreg,n_t,options,reduced_f,stochastic)
% IN:
%   - nreg: # regions in the network
%   - n_t: length of the time series
%   - options: options structure, containing at least the .inF and .inG
%   fields, as built using prepare_fullDCM.m
%   - reduced_f: a flag for fixing a subpart of the hemodynamic parameters
%   (oxygen extraction fraction at rest, vasodilatory signal feedback rate
%   and vessel stifness) to their prior value
%   - stochastic: a flag for stochastic DCM
% OUT:
%   - priors: the 'priors' structure that can be used to invert the DCM
%   using VBA_NLStateSpaceModel.m
%------------------------------------------------------------
% Copyright (C) 2012 Jean Daunizeau / License GNU GPL v2
%------------------------------------------------------------

dim.n_theta = options.inF.ind5(end);
if options.inG.homogeneous
    dim.n_phi = 2;
else
    dim.n_phi = 2*nreg;
end
dim.n = 5*nreg;

% initial conditions
priors.muX0 = zeros(5*nreg,1);
priors.SigmaX0 = 1e-4*eye(5*nreg);

% evolution parameters
priors.muTheta = zeros(dim.n_theta,1);
priors.SigmaTheta = 1e-3*eye(dim.n_theta);
priors.SigmaTheta(1:options.inF.indself,1:options.inF.indself) = ...
    1e-1*eye(options.inF.indself);
% priors.SigmaTheta(options.inF.indself,options.inF.indself) = 0;
if reduced_f
    % fix some HRF params to their default values
    priors.SigmaTheta(options.inF.ind1,options.inF.ind1) = 0;
    priors.SigmaTheta(options.inF.ind3,options.inF.ind3) = 0;
    priors.SigmaTheta(options.inF.ind5,options.inF.ind5) = 0;
end

% observation parameters
priors.muPhi = zeros(dim.n_phi,1);
priors.SigmaPhi = 1e-3*eye(dim.n_phi);

% state and measurement noise covariances
for t = 1:n_t
    dq = 1e2*ones(dim.n,1);
    dq(options.inF.n5) = 1;
    priors.iQx{t} = diag(dq);
    priors.iQy{t} = eye(nreg);
end
priors.AR = 0;

% precision hyperparameters
priors.a_sigma = 1e0;
priors.b_sigma = 1e0;
if ~stochastic
    priors.a_alpha = Inf;
    priors.b_alpha = 0;
else
%     TR = options.decim.*options.inF.deltat;
%     priors.b_alpha = TR/1e2;
    priors.a_alpha = 1e0;
    priors.b_alpha = 5e-2;
end




