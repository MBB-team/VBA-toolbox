function [Iphi,SigmaPhi,deltaMuPhi,suffStat] = VBA_Iphi_binomial(phi,y,posterior,suffStat,dim,u,options)
% Gauss-Newton update of the observation parameters, for binomial data

% check if called during initialization
if isequal(suffStat,VBA_getSuffStat(options))
    init = 1;
    if ~options.OnLine && options.verbose
        fprintf(1,'Deriving prior''s sufficient statistics ...')
        fprintf(1,'%6.2f %%',0)
    end
else
    init = 0;
end

if options.DisplayWin && ~init % Display progress
    if isequal(options.g_fname,@VBA_odeLim)
        STR = 'VB Gauss-Newton on observation/evolution parameters... ';
    else
        STR = 'VB Gauss-Newton on observation parameters... ';
    end
    set(options.display.hm(1),'string',STR);
    set(options.display.hm(2),'string','0%');
    drawnow
end

% Clear persistent variables if ODE mode
if isequal(options.g_fname,@VBA_odeLim) || ...
        isequal(options.g_fname,@VBA_smoothNLSS)
    clear VBA_odeLim
    clear VBA_smoothNLSS
end

%  Look-up which evolution parameter to update
indIn = options.params2update.phi;

% Preallocate intermediate variables
Q = options.priors.SigmaPhi(indIn,indIn);
iQ = VB_inv(Q,[]);
muPhi0 = options.priors.muPhi;
Phi = muPhi0;
Phi(indIn) = phi;
dphi0 = muPhi0-Phi;
ddydphi = 0;
dy = zeros(dim.p,dim.n_t);
vy = zeros(dim.p,dim.n_t);
gx = zeros(dim.p,dim.n_t);
d2gdx2 = zeros(dim.n_phi,dim.n_phi);
logL = 0;
if isequal(options.g_fname,@VBA_odeLim)
    muX = zeros(options.inG.old.dim.n,dim.n_t);
    SigmaX = cell(dim.n_t,1);
end
div = 0;

%--- Loop over time series ---%
for t=1:dim.n_t
    
    % evaluate observation function at current mode
    [gx(:,t),dG_dX,dG_dPhi,d2G_dXdPhi] = VBA_evalFun('g',posterior.muX(:,t),Phi,u(:,t),options,dim,t);
    
    % fix numerical instabilities
    gx(:,t) = checkGX_binomial(gx(:,t));
    
    % store states dynamics if ODE mode
    if isequal(options.g_fname,@VBA_odeLim)
        % get sufficient statistics of the hidden states from unused i/o in
        % VBA_evalFun.
        muX(:,t) = dG_dX;
        SigmaX{t} = d2G_dXdPhi'*posterior.SigmaPhi*d2G_dXdPhi;
    end
    
    % predicted variance over binomial data
    vy(:,t) = gx(:,t).*(1-gx(:,t));
    
    % remove irregular trials
    yin = find(~options.isYout(:,t));
    
    % accumulate log-likelihood
    logL = logL + y(yin,t)'*log(gx(yin,t)) + (1-y(yin,t))'*log(1-gx(yin,t));
    
    % prediction error
    dy(yin,t) = y(yin,t) - gx(yin,t);
    
    % gradient and Hessian
    ddydphi = ddydphi + dG_dPhi(:,yin)*(dy(yin,t)./vy(yin,t));
    tmp = y(yin,t)./gx(yin,t).^2 - (y(yin,t)-1)./(1-gx(yin,t)).^2;
    d2gdx2 = d2gdx2 + dG_dPhi(:,yin)*diag(tmp)*dG_dPhi(:,yin)';
    
    
    % Display progress
    if mod(t,dim.n_t./10) < 1
        if ~init && options.DisplayWin
            set(options.display.hm(2),'string',[num2str(floor(100*t/dim.n_t)),'%']);
            drawnow
        end
        if init && ~options.OnLine && options.verbose
            fprintf(1,repmat('\b',1,8))
            fprintf(1,'%6.2f %%',100*t/dim.n_t)
        end
    end
    
    % Accelerate divergent update
    if isweird({dy,dG_dPhi,dG_dX})
        div = 1;
        break
    end
    
end

% Display progress
if ~init && options.DisplayWin
    set(options.display.hm(2),'string','OK');
    drawnow
end
if init &&  ~options.OnLine  && options.verbose
    fprintf(1,repmat('\b',1,8))
    fprintf(' OK.')
    fprintf('\n')
end

% posterior covariance matrix
iSigmaPhi = iQ + d2gdx2(indIn,indIn);
SigmaPhi = VB_inv(iSigmaPhi,[]);

% mode
tmp = iQ*dphi0(indIn) + ddydphi(indIn);
deltaMuPhi = SigmaPhi*tmp;

% variational energy
Iphi = -0.5.*dphi0(indIn)'*iQ*dphi0(indIn) + logL;
if isweird({Iphi,SigmaPhi}) || div
    Iphi = -Inf;
end

% update sufficient statistics
suffStat.Sphi = 0.5*length(indIn)*log(2*pi*exp(1)) + 0.5*VBA_logDet(posterior.SigmaPhi,indIn);
suffStat.gx = gx;
suffStat.dy = dy;
suffStat.logL = logL;
suffStat.vy = vy;
suffStat.dphi = dphi0;
if isequal(options.g_fname,@VBA_odeLim)
    suffStat.muX = muX;
    suffStat.SigmaX = SigmaX;
end
suffStat.div = div;


