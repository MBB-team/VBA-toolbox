function [F] = VBA_FreeEnergy2(posterior,suffStat,options)
% computes free energy of the sDCM generative model (at equilibrium)
% function [F] = VBA_FreeEnergy2(posterior,suffStat,options)
% This function evaluates the free energy associated with the nonlinear
% state-space model inverted by the main inversion routine
% VBA_NLStateSpaceModel.m.
% IN:
%   - posterior: a structure containing the natural parameters of the
%   marginal posterior pdf of the unknown variables of the model
%   - suffStat: a structure containing pre-calculated (sufficient
%   statistics) quantities associated required for the computation of the
%   free energy (such as derivatives of the evolution/observation functions
%   evaluated at the current mode)
%   - options: a structure variable containing optional parameters (such as
%   the priors structure)
% OUT:
%   - F: the free energy under the Laplace approximation (with correction
%   terms due to precision hyperparameters).
%------------------------------------------------------------
% Copyright (C) 2012 Jean Daunizeau / License GNU GPL v2
%------------------------------------------------------------

if options.DisplayWin % Display progress
    try
        set(options.display.hm(1),'string',...
            'Calculating Free Energy... ');
        set(options.display.hm(2),'string','0%');
        drawnow
    end
end

priors = options.priors;
dim = options.dim;

% Get common free energy terms
SSE = (posterior.a_sigma./posterior.b_sigma)*suffStat.dy2;
dF = priors.a_sigma*log(priors.b_sigma) ...
    - posterior.a_sigma*log(posterior.b_sigma) ...
    + gammaln(posterior.a_sigma) - gammaln(priors.a_sigma) ...
    + posterior.a_sigma*(1-(priors.b_sigma/posterior.b_sigma));
ldQ = 0;
S = 0;
ntot = 0;

if dim.n > 0 && ~isinf(priors.a_alpha) && ~isequal(priors.b_alpha,0)
    SSE = SSE + (posterior.a_alpha./posterior.b_alpha)*suffStat.dx2;
    dF = dF + ...
        priors.a_alpha*log(priors.b_alpha) ...
        - posterior.a_alpha*log(posterior.b_alpha) ...
        + gammaln(posterior.a_alpha) - gammaln(priors.a_alpha) ...
        + posterior.a_alpha*(1-(priors.b_alpha/posterior.b_alpha));
    S = S + suffStat.SX;
end

for t=1:dim.n_t
    ldQ = ldQ + VBA_logDet(options.priors.iQy{t});
    ntot = ntot + length(find(diag(options.priors.iQy{t})~=0));
    if dim.n > 0  && ~isinf(priors.a_alpha) && ~isequal(priors.b_alpha,0)
        indIn = options.params2update.x{t};
        ldQ = ldQ + VBA_logDet(options.priors.iQx{t},indIn);
        ntot = ntot + length(indIn);
        S = S - 0.5*length(indIn);
    end
end


% observation parameters
if dim.n_phi > 0
    indIn = options.params2update.phi;
    if ~isempty(indIn)
        ntot = ntot + length(indIn);
        Q = priors.SigmaPhi(indIn,indIn);
        iQ = VB_inv(Q,[]);
        SSE = SSE + suffStat.dphi(indIn)'*iQ*suffStat.dphi(indIn);
        ldQ = ldQ - VBA_logDet(Q,[]);
        S = S + suffStat.Sphi - 0.5*length(indIn);
    end
end

% evolution parameters
if dim.n_theta > 0
    indIn = options.params2update.theta;
    if ~isempty(indIn)
        ntot = ntot + length(indIn);
        Q = priors.SigmaTheta(indIn,indIn);
        iQ = VB_inv(Q,[]);
        SSE = SSE + suffStat.dtheta(indIn)'*iQ*suffStat.dtheta(indIn);
        ldQ = ldQ - VBA_logDet(Q,[]);
        S = S + suffStat.Stheta - 0.5*length(indIn);
    end
end

% initial conditions
if dim.n > 0
    indIn = options.params2update.x0;
    if ~isempty(indIn)
        ntot = ntot + length(indIn);
        Q = priors.SigmaX0(indIn,indIn);
        iQ = VB_inv(Q,[]);
        SSE = SSE + suffStat.dx0(indIn)'*iQ*suffStat.dx0(indIn);
        ldQ = ldQ - VBA_logDet(Q,[]);
        S = S + suffStat.SX0 - 0.5*length(indIn);
    end
end


% Compose free energy
F = - 0.5*SSE - 0.5*ntot*log(2*pi) + 0.5*ldQ + S + dF;


if options.DisplayWin % Display progress
    try
        set(options.display.hm(2),'string','OK');
        drawnow
    end
end

