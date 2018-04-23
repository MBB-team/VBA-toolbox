function [options,dim] = getOptions4dcm(A,B,C,D,TR,microDT,n_t,homogeneous,reduced_f,lin)
% builds the options structure from basic DCM info
%
% function [options,dim] = getOptions4dcm(A,B,C,D,TR,microDT,n_t,homogeneous)
% IN:
%   - A: binary matrix indicating where the connections are
%   - B: cell-array of binary matrices of modulatory effects
%   - C: binary matrix of input-state coupling
%   - D: cell-array of binay matrices for gating effects
%   - TR: fMRI reptition time (data sampling period)
%   - microDT: micro-time resolution for ODE integration
%   - n_t: the number of time samples in the fMRI time series
%   - homogeneous: flag for indicating whether the observation parameters
%   of the Ballon model are identical across ROIs
%   - reduced_f: if 1, sets prior variance of 3 HRF parameters to 0 (fix
%   them to their prior mean)
%   - lin: a flag for using the linearized Balloon model
% OUT:
%   - options: optional structure for VB inversion of the specified model.

nreg = size(A,1);

try, homogeneous; catch, homogeneous = 1; end
try, reduced_f; catch, reduced_f = 1; end
try, lin; catch, lin = 0; end


[options]           = prepare_fullDCM(A,B,C,D,TR,microDT,homogeneous);
dim.n_theta         = options.inF.ind5(end);
dim.n_phi           = options.inG.ind2(end);
dim.n               = 5*nreg;
dim.p               = nreg;
dim.n_t             = n_t;



%--- Build priors for model inversion
thetaHRF = zeros(6,1);
priors.muX0 = kron(ones(nreg,1),[0;0;0;0;0]);
priors.SigmaX0 = 0e-1*eye(5*nreg);
priors.muTheta = 0e-1*ones(dim.n_theta,1);
priors.muTheta(options.inF.indself) = -0;
priors.muTheta(options.inF.ind1) = thetaHRF(1);
priors.muTheta(options.inF.ind3) = thetaHRF(3);
priors.muTheta(options.inF.ind5) = thetaHRF(6);
% priors.SigmaTheta = 1e-1*eye(dim.n_theta);
priors.SigmaTheta = 1e-2*eye(dim.n_theta);
priors.SigmaTheta(1:options.inF.indself,1:options.inF.indself) = 1e-1*eye(options.inF.indself);
priors.SigmaTheta(options.inF.indself,options.inF.indself) = 0;
if reduced_f
    % fix some HRF params to their default values
    priors.SigmaTheta(options.inF.ind1,options.inF.ind1) = 0;
    priors.SigmaTheta(options.inF.ind3,options.inF.ind3) = 0;
    priors.SigmaTheta(options.inF.ind5,options.inF.ind5) = 0;
end
priors.muPhi = 0*ones(dim.n_phi,1);
% priors.SigmaPhi = 1e-1*eye(dim.n_phi);
priors.SigmaPhi = 1e-2*eye(dim.n_phi);

% NB on hyperpriors:
%   - fix state noise precision using high scale param
%   - use non-informative priors on the residual precision, with high
%   expectation.
% This is because of the first iteration of the hidden states posterior
% update, which has to deviate from the its prior predictive density (as
% derived from the deterministic inversion).
% The following iterations will then work with a realistic (expected)
% residual precision, and adapt.

SC = 1e8;
priors.a_alpha = Inf;
priors.b_alpha = 0;
priors.a_sigma = 1e0;
priors.b_sigma = 1e0;


options.priors = priors;
options.updateHP = 1;
options.inF.linearized = lin;

options.sources = struct('type', 0, 'out', 1:dim.p);
