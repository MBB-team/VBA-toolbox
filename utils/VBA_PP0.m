function [PP] = VBA_PP0(posterior,out)
% derives the posterior probability P(param==0|y)
% function [PP] = VBA_PP0(posterior,out)
% The posterior probability P(param==0|y) is derived from a Bayesian model
% comparison between the original (full) model and a reduced model that
% assumes a priori P(param==0)=1. Note: this model comparison relies upon
% Savage-Dickey ratios [see VBA_SavageDickey.m].
% IN:
%   - posterior/out: output structures of VBA_NLStateSpaceModel.m
% OUT:
%   - PP: structure with fields corresponding to parameter sets, i.e.:
%       .phi: P(phi==0|y)
%       .theta: P(theta==0|y)
%       .X0: P(x0==0|y)

PP.phi = NaN(out.dim.n_phi,1);
PP.theta = NaN(out.dim.n_theta,1);
PP.X0 = NaN(out.dim.n,1);

% observation parameters
for i=1:out.dim.n_phi
    pr2 = out.options.priors;
    pr2.muPhi(i) = 0; % set prior mean to 0
    pr2.SigmaPhi(i,:) = 0; % set prior covariances to 0
    pr2.SigmaPhi(:,i) = 0; % set prior covariances to 0
    [F2] = VBA_SavageDickey(posterior,out.options.priors,out.F,out.dim,pr2);
    PP.phi(i) = VBA_sigmoid(F2-out.F);
end

% evolution parameters
for i=1:out.dim.n_theta
    pr2 = out.options.priors;
    pr2.muTheta(i) = 0; % set prior mean to 0
    pr2.SigmaTheta(i,:) = 0; % set prior covariances to 0
    pr2.SigmaTheta(:,i) = 0; % set prior covariances to 0
    [F2] = VBA_SavageDickey(posterior,out.options.priors,out.F,out.dim,pr2);
    PP.theta(i) = VBA_sigmoid(F2-out.F);
end

% initial conditions
for i=1:out.dim.n
    pr2 = out.options.priors;
    pr2.muX0(i) = 0; % set prior mean to 0
    pr2.SigmaX0(i,:) = 0; % set prior covariances to 0
    pr2.SigmaX0(:,i) = 0; % set prior covariances to 0
    [F2] = VBA_SavageDickey(posterior,out.options.priors,out.F,out.dim,pr2);
    PP.X0(i) = VBA_sigmoid(F2-out.F);
end

