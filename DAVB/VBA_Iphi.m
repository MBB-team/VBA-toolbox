function [Iphi,SigmaPhi,deltaMuPhi,suffStat] = VBA_Iphi(phi,y,posterior,suffStat,dim,u,options)
% Gauss-Newton update of the observation parameters
% !! When the observation function is @VBA_odeLim, this Gauss-Newton update
% actually implements a gradient ascent on the variational energy of the
% equivalent deterministic DCM.
%------------------------------------------------------------
% Copyright (C) 2012 Jean Daunizeau / License GNU GPL v2
%------------------------------------------------------------

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
if isequal(options.g_fname,@VBA_odeLim) || isequal(options.g_fname,@VBA_smoothNLSS)
    clear VBA_odeLim
    clear VBA_smoothNLSS
end

%  Look-up which evolution parameter to update
indIn = options.params2update.phi;

% Get precision parameters
sigmaHat = posterior.a_sigma./posterior.b_sigma;

% Preallocate intermediate variables
iQy = options.priors.iQy;
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
dy2 = 0;
Sphid2gdphi2 = 0;
kernel = zeros(dim.n_phi,dim.n_phi);
d2gdx2 = 0;
if isequal(options.g_fname,@VBA_odeLim)
    muX = zeros(options.inG.old.dim.n,dim.n_t);
    SigmaX = cell(dim.n_t,1);
end
div = 0;

%--- Loop over time series ---%
for t=1:dim.n_t
    
    % evaluate observation function at current mode
    [gx(:,t),dG_dX,dG_dPhi,d2G_dXdPhi] = VBA_evalFun('g',posterior.muX(:,t),Phi,u(:,t),options,dim);
    
    if isequal(options.g_fname,@VBA_odeLim)
        % get sufficient statistics of the hidden states from unused i/o in
        % VBA_evalFun.
        muX(:,t) = dG_dX;
        SigmaX{t} = d2G_dXdPhi'*posterior.SigmaPhi*d2G_dXdPhi;
    end
    
    % mean-field terms
    Sphid2gdphi2 = Sphid2gdphi2 + trace(dG_dPhi*iQy{t}*dG_dPhi'*posterior.SigmaPhi);
    
    % posterior covariance matrix terms
    d2gdx2 = d2gdx2 + dG_dPhi*iQy{t}*dG_dPhi';
    
    % error terms
    dy(:,t) = y(:,t) - gx(:,t);
    dy2 = dy2 + dy(:,t)'*iQy{t}*dy(:,t);
    ddydphi = ddydphi + dG_dPhi*iQy{t}*dy(:,t);
    if dim.n > 0 && ~options.ignoreMF
        A1g = reshape(permute(d2G_dXdPhi,[1,3,2]),dim.p*dim.n,dim.n_phi)';
        A2g = A1g*kron(iQy{t},posterior.SigmaX.current{t});
        kernel = kernel + A2g*A1g';
        ddydphi = ddydphi + A2g*vec(dG_dX);
    end
    
    % Predictive density (data space)
    V = dG_dPhi'*posterior.SigmaPhi*dG_dPhi + (1./sigmaHat).*VB_inv(iQy{t},[]);
    if dim.n > 0
        V = V + dG_dX'*posterior.SigmaX.current{t}*dG_dX;
    end
    vy(:,t) = diag(V);
    
    
    % Display progress
    if mod(t,dim.n_t./10) < 1
        if ~init && options.DisplayWin
            set(options.display.hm(2),...
                'string',[num2str(floor(100*t/dim.n_t)),'%']);
            drawnow
        end
        if init && ~options.OnLine && options.verbose
            fprintf(1,repmat('\b',1,8))
            fprintf(1,'%6.2f %%',100*t/dim.n_t)
        end
    end
    
    % Accelerate divergent update
    if isweird({dy2,dG_dPhi,dG_dX})
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
    fprintf(1,' OK.')
    fprintf(1,'\n')
end


% posterior covariance matrix
iSigmaPhi = iQ + sigmaHat.*( d2gdx2(indIn,indIn) + kernel(indIn,indIn) );
SigmaPhi = VB_inv(iSigmaPhi,[]);

% mode
tmp = iQ*dphi0(indIn) + sigmaHat.*ddydphi(indIn);
deltaMuPhi = SigmaPhi*tmp;

% variational energy
Iphi = -0.5.*dphi0(indIn)'*iQ*dphi0(indIn) -0.5*sigmaHat.*dy2;
if ~options.ignoreMF
    Iphi = Iphi -0.5*sigmaHat.*SXd2gdx2;
end
if isweird({Iphi,SigmaPhi}) || div
    Iphi = -Inf;
end

% update sufficient statistics
suffStat.Sphi = 0.5*length(indIn)*log(2*pi*exp(1)) + 0.5*VBA_logDet(posterior.SigmaPhi,indIn);
suffStat.Sphid2gdphi2 = Sphid2gdphi2;
suffStat.Sphid2gdphidx = trace(kernel*posterior.SigmaPhi);
suffStat.gx = gx;
suffStat.dy = dy;
suffStat.dy2 = dy2;
suffStat.vy = vy;
suffStat.dphi = dphi0;
if isequal(options.g_fname,@VBA_odeLim)
    suffStat.muX = muX;
    suffStat.SigmaX = SigmaX;
end
suffStat.div = div;



