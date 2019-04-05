function [Iphi,SigmaPhi,deltaMuPhi,suffStat] = VBA_Iphi_split(phi,y,posterior,suffStat,dim,u,options)
% Gauss-Newton update of the observation parameters
% !! When the observation function is @VBA_odeLim, this Gauss-Newton update
% actually implements a gradient ascent on the variational energy of the
% equivalent deterministic DCM.

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
if isequal(options.g_fname,@VBA_odeLim)
    clear VBA_odeLim
    muX = zeros(options.inG.old.dim.n,dim.n_t);
    SigmaX = cell(dim.n_t,1);
end

%  Look-up which evolution parameter to update
indIn = options.params2update.phi;

% Get precision parameters
sigmaHat = posterior.a_sigma./posterior.b_sigma;

% Preallocate intermediate variables
iQy = options.priors.iQy;
Q = options.priors.SigmaPhi(indIn,indIn);
iQ = VBA_inv(Q);
muPhi0 = options.priors.muPhi;
Phi = muPhi0;
Phi(indIn) = phi;
ddydphi = 0;
dy = zeros(dim.p,dim.n_t);
vy = zeros(dim.p,dim.n_t);
gx = zeros(dim.p,dim.n_t);
dy2 = 0;
d2gdx2 = 0;
div = 0;

% intermediary variables: MoG split
sqrtS = VBA_sqrtm (posterior.SigmaPhi(indIn,indIn));
split = options.split;
nd = size(split.m,2);
Mu0 = muPhi0;
Sigma0 = options.priors.SigmaPhi;

%--- Loop over time series ---%
for t=1:dim.n_t
    
    if isequal(options.g_fname,@VBA_odeLim)
        SigmaX{t} = 0;
        muXt = zeros(options.inG.old.dim.n,nd);
    end
    Vy = cell(nd,1);
    Vx = cell(nd,1);
    gxt = zeros(dim.p,nd);
    
    % loop over n_phi-draws
    for i=1:nd
        
        Mu0(indIn) = sqrtS*split.m(:,i) + Phi;
        Sigma0(indIn,indIn) = sqrtS*diag(split.s(:,i))*sqrtS';
        
        % evaluate observation function at current mode
        [gxt(:,i),dG_dX,dG_dPhi] = VBA_evalFun('g',posterior.muX(:,t),Mu0,u(:,t),options,dim,t);
        
        % mean-field terms
        if isequal(options.g_fname,@VBA_odeLim)
            % get sufficient statistics of the hidden states from unused i/o in
            % VBA_evalFun.
            muXt(:,i) = dG_dX.xt;
            Vx{i} = dG_dX.dx'*Sigma0*dG_dX.dx;
        end
        
        % posterior covariance matrix terms
        d2gdx2 = d2gdx2 + split.w(i).*dG_dPhi*iQy{t}*dG_dPhi';
        
        % error terms
        dyti = y(:,t) - gxt(:,i);
        dy(:,t) = dy(:,t)+ split.w(i).*dyti;
        dy2 = dy2 + split.w(i).*dyti'*iQy{t}*dyti;
        ddydphi = ddydphi + split.w(i).*dG_dPhi*iQy{t}*dyti;
        
        % Predictive density (data space)
        Vy{i} = dG_dPhi'*Sigma0*dG_dPhi + (1./sigmaHat).*VBA_inv(iQy{t});
        if dim.n > 0
            Vy{i} = Vy{i} + dG_dX'*posterior.SigmaX.current{t}*dG_dX;
        end
        
    end
    
    % 4- match moments of the MoG
    gx(:,t) = sum(gxt*diag(split.w(:)),2);
    if isequal(options.g_fname,@VBA_odeLim)
        muX(:,t) = sum(muXt*diag(split.w(:)),2);
    end
    V = 0;
    for i=1:nd
        dg = gxt(:,i) - gx(:,t);
        V = V + split.w(i).*(dg*dg'+Vy{i});
        if isequal(options.g_fname,@VBA_odeLim)
            dxp = muXt(:,i) - muX(:,t);
            SigmaX{t} = SigmaX{t} + split.w(i).*(dxp*dxp'+Vx{i});
        end
    end
    vy(:,t) = diag(V);
    
    
    % Display progress
    if isequal(mod(t,dim.n_t./10),0)
        if ~init && options.DisplayWin
            set(options.display.hm(2),'string',[num2str(100*t/dim.n_t),'%']);
            drawnow
        end
        if init && ~options.OnLine && options.verbose
            fprintf(1,repmat('\b',1,8))
            fprintf(1,'%6.2f %%',100*t/dim.n_t)
        end
    end
    
    % Accelerate divergent update
    if VBA_isWeird ({dy2, dG_dPhi, dG_dX})
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
iSigmaPhi = iQ + sigmaHat.*d2gdx2(indIn,indIn);
SigmaPhi = VBA_inv(iSigmaPhi);%./split.s(1);

% mode
Mu0 = muPhi0;
dphi0 = 0;
for i=1:nd
    Mu0(indIn) = sqrtS*split.m(:,i) + Phi;
    dphi0 = dphi0 + split.w(i).*(muPhi0-Mu0);
end
tmp = iQ*dphi0(indIn) + sigmaHat.*ddydphi(indIn);
deltaMuPhi = SigmaPhi*tmp;%split.s(1).*SigmaPhi*tmp;

% variational energy
Iphi = -0.5.*dphi0(indIn)'*iQ*dphi0(indIn) -0.5*sigmaHat.*dy2;
if VBA_isWeird ({Iphi, SigmaPhi}) || div
    Iphi = -Inf;
end

% update sufficient statistics
suffStat.gx = gx;
suffStat.dy = dy;
suffStat.dy2 = dy2;
suffStat.vy = vy;
suffStat.dphi = muPhi0 - Phi;
if isequal(options.g_fname,@VBA_odeLim)
    suffStat.muX = muX;
    suffStat.SigmaX = SigmaX;
end



