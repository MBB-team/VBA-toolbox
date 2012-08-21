function [F] = VBA_FreeEnergy(posterior,suffStat,options)
% computes free energy of the sDCM generative model
% function [F] = VBA_FreeEnergy(posterior,suffStat,options)
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
%   - F: the free energy under the local Laplace approximation



priors = options.priors;
dim = options.dim;

if options.binomial
    F = suffStat.logL;
    % observation parameters
    if dim.n_phi > 0
        indIn = options.params2update.phi;
        if ~isempty(indIn)
            Q = priors.SigmaPhi(indIn,indIn);
            iQ = VB_inv(Q,[]);
            F = F ...
                - 0.5*suffStat.dphi(indIn)'*iQ*suffStat.dphi(indIn) ...
                - 0.5*length(indIn)*log(2*pi) ...
                - 0.5*VBA_logDet(Q,[]) ...
                + suffStat.Sphi - 0.5*length(indIn);
        end
    end
    return
end


if options.Laplace
    [F] = VBA_FreeEnergy3(posterior,suffStat,options);
%     [dF] = F-VBA_FreeEnergy2(posterior,suffStat,options)
    return
end

if options.DisplayWin % Display progress
    try
        set(options.display.hm(1),'string',...
            'Calculating Free Energy... ');
        set(options.display.hm(2),'string','0%');
        drawnow
    end
end


% Get precision parameters
sigmaHat = posterior.a_sigma./posterior.b_sigma;
logSigmaHat = VBA_psi(posterior.a_sigma) - log(posterior.b_sigma);
if dim.n > 0
    alphaHat = posterior.a_alpha./posterior.b_alpha;
    logAlphaHat = VBA_psi(posterior.a_alpha) - log(posterior.b_alpha);
end

ldQ = 0;
nnt = 0;
for t=1:dim.n_t
    ldQ = ldQ + VBA_logDet(options.priors.iQy{t});
    if dim.n > 0
        ldQ = ldQ + ...
            VBA_logDet(options.priors.iQx{t},options.params2update.x{t});
        nnt = nnt + length(options.params2update.x{t});
    end
end


% Common terms in the free energy
F = - 0.5*sigmaHat*suffStat.dy2 ...
    - 0.5*dim.p*dim.n_t*(log(2*pi) - logSigmaHat) ...
    + priors.a_sigma*log(priors.b_sigma) - gammaln(priors.a_sigma) ...
    + (priors.a_sigma-1).*logSigmaHat - sigmaHat.*priors.b_sigma ...
    + suffStat.Ssigma ...
    + 0.5*ldQ;

% observation parameters
if dim.n_phi > 0
    indIn = options.params2update.phi;
    if ~isempty(indIn)
        Q = priors.SigmaPhi(indIn,indIn);
        Sphi = posterior.SigmaPhi(indIn,indIn);
        iQ = VB_inv(Q,[]);
        F = F ...
            - 0.5*sigmaHat.*suffStat.Sphid2gdphi2 ...
            - 0.5*suffStat.dphi(indIn)'*iQ*suffStat.dphi(indIn) ...
            - 0.5.*trace(iQ*Sphi) ...
            - 0.5*length(indIn)*log(2*pi) ...
            - 0.5*VBA_logDet(Q) ...
            + suffStat.Sphi;
    end
    if dim.n > 0
        F = F ...
            - 0.5*sigmaHat.*suffStat.Sphid2gdphidx;
    end
end

% evolution parameters
if dim.n_theta > 0
    indIn = options.params2update.theta;
    if ~isempty(indIn)
        Q = priors.SigmaTheta(indIn,indIn);
        Stheta = posterior.SigmaTheta(indIn,indIn);
        iQ = VB_inv(Q,[]);
        F = F ...
            - 0.5*alphaHat.*suffStat.Sthetad2fdtheta2 ...
            - 0.5*suffStat.dtheta(indIn)'*iQ*suffStat.dtheta(indIn) ...
            - 0.5.*trace(iQ*Stheta) ...
            - 0.5*length(indIn)*log(2*pi) ...
            - 0.5*VBA_logDet(Q) ...
            + suffStat.Stheta;
    end
    if dim.n > 0
        F = F ...
            - 0.5*alphaHat.*suffStat.Sthetad2fdthetadx;
    end
end

% hidden states
if dim.n > 0
    F = F ...
        - 0.5*alphaHat*suffStat.dx2 ...
        - 0.5*nnt*(log(2*pi) - logAlphaHat) ...
        + priors.a_alpha*log(priors.b_alpha) - gammaln(priors.a_alpha) ...
        + (priors.a_alpha-1).*logAlphaHat - alphaHat.*priors.b_alpha ...
        + suffStat.Salpha ...
        - 0.5*alphaHat.*suffStat.SXd2fdx2 ...
        - 0.5*alphaHat.*suffStat.SXtdfdx ...
        - 0.5*sigmaHat.*suffStat.SXd2gdx2 ...
        - 0.5*alphaHat.*suffStat.trSx ...
        + suffStat.SX;
    if options.updateX0
        indIn = options.params2update.x0;
        Q = priors.SigmaX0(indIn,indIn);
        SX0 = posterior.SigmaX0(indIn,indIn);
        iQ = VB_inv(Q,[]);
        F = F ...
            - 0.5*suffStat.dx0(indIn)'*iQ*suffStat.dx0(indIn) ...
            - 0.5.*trace(iQ*SX0) ...
            - 0.5*length(indIn)*log(2*pi) ...
            - 0.5*VBA_logDet(Q) ...
            + suffStat.SX0;
    end
end

if options.DisplayWin % Display progress
    try
        set(options.display.hm(2),'string','OK');
        drawnow
    end
end
