function suffStat = VBA_Hpost(posterior,suffStat,options)
% Entropy calculus

if options.dim.n > 0
    % initial conditions
    indIn = options.params2update.x0;
    suffStat.SX0 = 0.5*length(indIn)*log(2*pi*exp(1)) + 0.5*VBA_logDet(posterior.SigmaX0,indIn);
    % hidden states
    indIn = options.params2update.x;
    SX = 0.5*length(indIn{1})*log(2*pi*exp(1)) + 0.5*VBA_logDet(posterior.SigmaX.current{1},indIn{1});
    for t = 2:options.dim.n_t-1
        jointCov = ...
            [   posterior.SigmaX.current{t+1}   posterior.SigmaX.inter{t}'
                posterior.SigmaX.inter{t}       posterior.SigmaX.current{t} ];
        indjc = [indIn{t+1}(:);indIn{t}(:)+options.dim.n];
        ldj = VBA_logDet(jointCov(indjc,indjc));
        ldm = VBA_logDet(posterior.SigmaX.current{t}(indIn{t},indIn{t}));
        SX = SX + 0.5*( ldj - ldm ) + 0.5*length(indIn{t})*log(2*pi*exp(1));
    end
    % precision hyperparameter
    if ~isinf(options.priors.a_alpha) && ~isequal(options.priors.b_alpha,0)
        suffStat.Salpha = gammaln(posterior.a_alpha) ...
            - log(posterior.b_alpha) ...
            + (1-posterior.a_alpha).*psi(posterior.a_alpha) ...
            + posterior.a_alpha;
    end
end

if options.dim.n_theta > 0
    indIn = options.params2update.theta;
    suffStat.Stheta = 0.5*length(indIn)*log(2*pi*exp(1)) + 0.5*VBA_logDet(posterior.SigmaTheta,indIn);
end

if options.dim.n_phi > 0
    indIn = options.params2update.phi;
    suffStat.Sphi = 0.5*length(indIn)*log(2*pi*exp(1)) + 0.5*VBA_logDet(posterior.SigmaPhi,indIn);
end

if ~options.binomial
    suffStat.Ssigma = gammaln(posterior.a_sigma) ...
        - log(posterior.b_sigma) ...
        + (1-posterior.a_sigma).*psi(posterior.a_sigma) ...
        + posterior.a_sigma;
end

