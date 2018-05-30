function [F] = VBA_FreeEnergy_UNL(posterior,suffStat,options)
% computes free energy of generative models with unnormalized likelihoods
% function [F] = VBA_FreeEnergy_UNL(posterior,suffStat,options)
% IN:
%   - posterior: a structure containing the natural parameters of the
%   marginal posterior pdf of the unknown model variables
%   - suffStat: a structure containing pre-calculated (sufficient
%   statistics) quantities associated required for the computation of the
%   free energy (such as derivatives of the evolution/observation functions
%   evaluated at the current mode)
%   - options: a structure variable containing optional parameters (such as
%   the priors structure)
% OUT:
%   - F: the free energy under the local Laplace approximation

if options.DisplayWin % Display progress
    try
        set(options.display.hm(1),'string','Calculating Free Energy... ');
        set(options.display.hm(2),'string','0%');
        drawnow
    end
end

priors = options.priors;
dim = options.dim;

% Entropy calculus
indIn = options.params2update.phi;
suffStat.Sphi = 0.5*length(indIn)*log(2*pi*exp(1)) + 0.5*VBA_logDet(posterior.SigmaPhi,indIn);
F = sum(suffStat.logL);

% observation parameters
if dim.n_phi > 0
    indIn = options.params2update.phi;
    if ~isempty(indIn)
        ntot = length(indIn);
        Q = priors.SigmaPhi(indIn,indIn);
        iQ = VBA_inv(Q);
        SSE = suffStat.dphi(indIn)'*iQ*suffStat.dphi(indIn);
        ldQ = - VBA_logDet(Q,[]);
        S = suffStat.Sphi - 0.5*length(indIn);
    end
    F = F - 0.5*SSE - 0.5*ntot*log(2*pi) + 0.5*ldQ + S;
end

E = posterior.a_sigma./posterior.b_sigma;
V = posterior.a_sigma./posterior.b_sigma^2;
E0 = priors.a_sigma./priors.b_sigma;
V0 = priors.a_sigma./priors.b_sigma^2;
F = F -VBA_KL(E,V,E0,V0,'Gamma');

if options.DisplayWin % Display progress
    try
        set(options.display.hm(2),'string','OK');
        drawnow
    end
end

