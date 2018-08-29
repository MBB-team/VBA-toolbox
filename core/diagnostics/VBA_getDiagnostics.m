function [diagnostics,out] = VBA_getDiagnostics(posterior,out)
% derives post-hoc diagnostics of VBA's model inversion
% function [diagnostics,out] = VBA_getDiagnostics(posterior,out)
% Post-hoc diagnostics include: goodness-of-fit metrics, Volterra kernels,
% null model evidence, posterior entropies, etc...
% IN:
%   - posterior,out: VBA's output structures
% OUT:
%   - diagnostics: a structure containing the following fields:
%       .pgx: mean of the prior predictive density
%       .pvy: variance of the prior predictive density
%       .kernels: Volterra kernels (if options.kernelSize>0)
%       .efficiency: posterior entropies
%       .DKL: prior/posterior Kullback-Leibler divergences
%       .LLH0: log-evidence of the null model
%       .MT_x: microtime hidden states time series
%       .MT_gx: microtime predicted data time series
%       .microTime: microtime sampling grid
%       .sampleInd: sub-indexing of the microtime sampling grid
%       .dy: data residuals structure (e.g., autocorrelation...)
%       .dx: states innovations structure
%       .C: posterior correlation matrix

u = out.u;
y = out.y;

% get goodness-tof-fit metrics
try; out.fit; catch; out.fit = VBA_getFit(posterior,out); end

% derive Volterra kernels
if out.dim.n_t>1 && out.options.kernelSize>0
    try
        kernels = VBA_getVolterraKernels(posterior,out);
    catch
        VBA_disp('  *** could not derive kernels!\n',out.options);
        kernels = [];
    end
else
    kernels = [];
end

% get null model (H0) evidence
[LLH0] = VBA_LMEH0(y,out.options);

% Entropies and KL divergences
gsi = find([out.options.sources.type]==0) ;
if isempty(gsi)
    efficiency.sigma = NaN;
    DKL.sigma = NaN;
else
    for iG=1:numel(gsi)
        efficiency.sigma(iG) = -out.suffStat.Ssigma(iG);
        m0 = out.options.priors.a_sigma(iG)/out.options.priors.b_sigma(iG);
        v0 = out.options.priors.a_sigma(iG)/out.options.priors.b_sigma(iG)^2;
        m = posterior.a_sigma(iG)/posterior.b_sigma(iG);
        v = posterior.a_sigma(iG)/posterior.b_sigma(iG)^2;
        DKL.sigma(iG) = VBA_KL(m0,v0,m,v,'Gamma');
    end
end
if out.dim.n > 0 % hidden states and initial conditions
    efficiency.X = -out.suffStat.SX;
    efficiency.X0 = -out.suffStat.SX0;
    if isinf(out.options.priors.a_alpha) && isequal(out.options.priors.b_alpha,0)
        efficiency.alpha = NaN;
        DKL.alpha = NaN;
    else
        efficiency.alpha = -out.suffStat.Salpha;
        m0 = out.options.priors.a_alpha./out.options.priors.b_alpha;
        v0 = out.options.priors.a_alpha./out.options.priors.b_alpha^2;
        m = posterior.a_alpha(end)./posterior.b_alpha(end);
        v = posterior.a_alpha(end)./posterior.b_alpha(end)^2;
        DKL.alpha = VBA_KL(m0,v0,m,v,'Gamma');
    end
    try
        DKL.X = 0;
        for t=1:out.dim.n_t
            IN = out.options.params2update.x{t};
            m0 = out.options.priors.muX(IN,t);
            v0 = out.options.priors.SigmaX.current{t}(IN,IN);
            m = posterior.muX(IN,t);
            v = posterior.SigmaX.current{t}(IN,IN);
            DKL.X = DKL.X + VBA_KL(m0,v0,m,v,'Normal');
        end
    catch
        DKL.X = NaN;
    end
    IN = out.options.params2update.x0;
    m0 = out.options.priors.muX0(IN);
    v0 = out.options.priors.SigmaX0(IN,IN);
    m = posterior.muX0(IN);
    v = posterior.SigmaX0(IN,IN);
    DKL.X0 = VBA_KL(m0,v0,m,v,'Normal');
else
    efficiency.X = NaN;
    efficiency.X0 = NaN;
    efficiency.alpha = NaN;
    DKL.X = NaN;
    DKL.X0 = NaN;
    DKL.alpha = NaN;
end
if out.dim.n_phi > 0 % observation parameters
    efficiency.Phi = -out.suffStat.Sphi;
    IN = out.options.params2update.phi;
    m0 = out.options.priors.muPhi(IN);
    v0 = out.options.priors.SigmaPhi(IN,IN);
    if ~out.options.OnLine
        m = posterior.muPhi(IN);
        v = posterior.SigmaPhi(IN,IN);
    else
        m = posterior.muPhi(IN,end);
        v = posterior.SigmaPhi{end}(IN,IN);
    end
    DKL.Phi = VBA_KL(m0,v0,m,v,'Normal');
else
    efficiency.Phi = NaN;
    DKL.Phi = NaN;
end
if out.dim.n_theta > 0 % evolution parameters
    efficiency.Theta = -out.suffStat.Stheta;
    IN = out.options.params2update.theta;
    m0 = out.options.priors.muTheta(IN);
    v0 = out.options.priors.SigmaTheta(IN,IN);
    if ~out.options.OnLine
        m = posterior.muTheta(IN);
        v = posterior.SigmaTheta(IN,IN);
    else
        m = posterior.muTheta(IN,end);
        v = posterior.SigmaTheta{end}(IN,IN);
    end
    DKL.Theta = VBA_KL(m0,v0,m,v,'Normal');
else
    efficiency.Theta = NaN;
    DKL.Theta = NaN;
end

% get prior predictive density
try
    [muy,Vy] = VBA_getLaplace(u,out.options.f_fname,out.options.g_fname,out.dim,out.options,0,'diag');
catch
    muy =[];
    Vy = [];
end

% get micro-time posterior hidden-states estimates
try
    [MT_x,MT_gx,microTime,sampleInd] = VBA_microTime(posterior,u,out);
catch
    MT_x = [];
    MT_gx = [];
    microTime = 1:out.dim.n_t;
    sampleInd = 1:out.dim.n_t;
end

% get residuals structure: data noise
for iS = 1:numel(out.options.sources)
    % source wise
    ySource = out.options.sources(iS).out ;
    res = out.suffStat.dy(ySource,:);
    if out.options.sources(iS).type==0
        gi = find(iS==find([out.options.sources(:).type]==0));
        iQyt = out.options.priors.iQy(:,gi);
        res = getWeightedResiduals(res,iQyt);
    end
    % remove skipped
    res(out.options.isYout(ySource,:)==1) = NaN ;
    dy(iS).dy = res(:);
    dy(iS).R = VBA_spm_autocorr(res);
    dy(iS).m = VBA_nanmean(dy(iS).dy);
    dy(iS).v = VBA_nanvar(dy(iS).dy);
    [dy(iS).ny,dy(iS).nx] = hist(dy(iS).dy,10);
    dy(iS).ny = dy(iS).ny./sum(dy(iS).ny);
    d = diff(dy(iS).nx);
    d = abs(d(1));
    dy(iS).d = d;
    spgy = sum(exp(-0.5.*(dy(iS).m-dy(iS).nx).^2./dy(iS).v));
    dy(iS).grid = dy(iS).nx(1):d*1e-2:dy(iS).nx(end);
    dy(iS).pg = exp(-0.5.*(dy(iS).m-dy(iS).grid).^2./dy(iS).v);
    dy(iS).pg = dy(iS).pg./spgy;
    if  out.options.sources(iS).type==0
        igs = sum([out.options.sources(1:iS).type]==0);
        shat = posterior.a_sigma(igs)./posterior.b_sigma(igs);
        spgy = sum(exp(-0.5.*shat.*dy(iS).nx.^2));
        dy(iS).pg2 = exp(-0.5.*shat.*dy(iS).grid.^2);
        dy(iS).pg2 = dy(iS).pg2./spgy;
    end
end

% get residuals structure: state noise
if ~isempty(out.suffStat.dx)
    wdx = getWeightedResiduals(out.suffStat.dx,out.options.priors.iQx);
    dx.dx = wdx(:);
    dx.m = mean(dx.dx);
    dx.v = var(dx.dx);
    [dx.ny,dx.nx] = hist(dx.dx,10);
    dx.ny = dx.ny./sum(dx.ny);
    d = diff(dx.nx);
    d = abs(d(1));
    dx.d = d;
    spgy = sum(exp(-0.5.*(dx.m-dx.nx).^2./dx.v));
    dx.grid = dx.nx(1):d*1e-2:dx.nx(end);
    dx.pg = exp(-0.5.*(dx.m-dx.grid).^2./dx.v);
    dx.pg = dx.pg./spgy;
    ahat = posterior.a_alpha(end)./posterior.b_alpha(end);
    spgy = sum(exp(-0.5.*ahat.*dx.nx.^2));
    dx.pg2 = exp(-0.5.*ahat.*dx.grid.^2);
    dx.pg2 = dx.pg2./spgy;
else
    dx.dx = [];
end

% get parameters posterior correlation matrix
if out.dim.n > 0 && isinf(out.options.priors.a_alpha) && isequal(out.options.priors.b_alpha,0)
    S = out.suffStat.ODE_posterior.SigmaPhi;
else
    S = NaN*zeros(out.dim.n+out.dim.n_theta+out.dim.n_phi);
    ind = 0;
    if out.dim.n_phi > 0
        if iscell(posterior.SigmaPhi) % online version
            SP = posterior.SigmaPhi{end};
        else
            SP = posterior.SigmaPhi;
        end
        S(1:out.dim.n_phi,1:out.dim.n_phi) = SP;
        ind = out.dim.n_phi;
    end
    if out.dim.n_theta > 0
        if iscell(posterior.SigmaTheta) % online version
            SP = posterior.SigmaTheta{end};
        else
            SP = posterior.SigmaTheta;
        end
        S(ind+1:ind+out.dim.n_theta,ind+1:ind+out.dim.n_theta) = SP;
        ind = ind + out.dim.n_theta;
    end
    if out.dim.n > 0 && out.options.updateX0
        if iscell(posterior.SigmaX0) % online version
            SP = posterior.SigmaX0{end};
        else
            SP = posterior.SigmaX0;
        end
        S(ind+1:ind+out.dim.n,ind+1:ind+out.dim.n) = SP;
    end
end
C = VBA_cov2corr(S);
C = C + diag(NaN.*diag(C));
tick = [0];
ltick = [];
ticklabel = cell(0,0);
if out.dim.n_phi > 0
    ltick = [ltick,tick(end)+out.dim.n_phi/2];
    tick = [tick,out.dim.n_phi];
    ticklabel{end+1} = 'phi';
end
if out.dim.n_theta > 0
    ltick = [ltick,tick(end)+out.dim.n_theta/2];
    tick = [tick,tick(end)+out.dim.n_theta];
    ticklabel{end+1} = 'theta';
end
if out.dim.n > 0 && out.options.updateX0
    ltick = [ltick,tick(end)+out.dim.n/2];
    tick = [tick,tick(end)+out.dim.n];
    ticklabel{end+1} = 'x0';
end
tick = tick +0.5;
tick = tick(2:end-1);
ltick = ltick + 0.5;


% wrap up
diagnostics.pgx = reshape(muy,out.dim.p,[]);
diagnostics.pvy = reshape(Vy,out.dim.p,[]);
diagnostics.kernels = kernels;
diagnostics.efficiency = efficiency;
diagnostics.DKL = DKL;
diagnostics.LLH0 = LLH0;
diagnostics.MT_x = MT_x;
diagnostics.MT_gx = MT_gx;
diagnostics.microTime = microTime;
diagnostics.sampleInd = sampleInd;
diagnostics.dy = dy;
diagnostics.dx = dx;
diagnostics.ltick = ltick;
diagnostics.tick = tick;
diagnostics.ticklabel = ticklabel;
diagnostics.C = C;
out.diagnostics = diagnostics;


function wdx = getWeightedResiduals(dx,iQx)
% weigths residuals according to (state/data) precision matrix
wdx = zeros(size(dx));
for t = 1:size(dx,2)
    sqrtiQ = VBA_sqrtm (iQx{t});
    wdx(:,t) = sqrtiQ*dx(:,t);
end



