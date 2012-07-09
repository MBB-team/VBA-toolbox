function [suffStat,posterior] = VBA_check_errors(y,u,options)
% Function merging VBA_getsuffstat and the content of VBA_Initialize
% regarding deterministic DCM
% ---------------------------------------------------------------------
% ------- VBA_getsuffstats
% [suffStat] = VBA_getSuffStat(options,suffStat,flag)
% fills in default sufficient statistics structure
% function [suffStat] = VBA_getSuffStat(options)
% IN:
%   - options: the options structure (only the substructure .dim is used
%   here)
%   - suffStat: {[]}
%   - flag: {0}, 1 if online version


flag = 0;
suffStat = [];

try
    flag;
catch
    flag = 0;
end

try
    suffStat.gx;
catch
    suffStat.gx = zeros(options.dim.p,options.dim.n_t);
end
try
    suffStat.vy;
catch
    suffStat.vy = zeros(options.dim.p,options.dim.n_t);
end
try
    suffStat.muX;
catch
    suffStat.muX = zeros(options.dim.n,options.dim.n_t);
end
try
    suffStat.dy;
catch
    suffStat.dy = zeros(options.dim.p,options.dim.n_t);
end
try
    suffStat.dx;
catch
    suffStat.dx = zeros(options.dim.n,options.dim.n_t);
end
try
    suffStat.vdx;
catch
    suffStat.vdx = zeros(options.dim.n,options.dim.n_t);
end
try
    suffStat.dx0;
catch
    suffStat.dx0 = zeros(options.dim.n,1);
end
try
    suffStat.dy2;
catch
    suffStat.dy2 = 0;
end
try
    suffStat.dx2;
catch
    suffStat.dx2 = 0;
end
try
    suffStat.dtheta;
catch
    if ~flag
        suffStat.dtheta = zeros(options.dim.n_theta,1);
    else
        suffStat.dtheta = zeros(options.dim.n_theta,options.dim.n_t);
    end
end
try
    suffStat.dphi;
catch
    if ~flag
        suffStat.dphi = zeros(options.dim.n_phi,1);
    else
        suffStat.dphi = zeros(options.dim.n_phi,options.dim.n_t);
    end
end
try
    suffStat.Ssigma;
catch
    suffStat.Ssigma = 0;
end
try
    suffStat.Sphid2gdphi2;
catch
    suffStat.Sphid2gdphi2 = 0;
end
try
    suffStat.Sphi;
catch
    suffStat.Sphi = 0;
end
try
    suffStat.Sphid2gdphidx;
catch
    suffStat.Sphid2gdphidx = 0;
end
try
    suffStat.Sthetad2fdtheta2;
catch
    suffStat.Sthetad2fdtheta2 = 0;
end
try
    suffStat.Stheta;
catch
    suffStat.Stheta = 0;
end
try
    suffStat.Sthetad2fdthetadx;
catch
    suffStat.Sthetad2fdthetadx = 0;
end
try
    suffStat.Salpha;
catch
    suffStat.Salpha = 0;
end
try
    suffStat.SXd2fdx2;
catch
    suffStat.SXd2fdx2 = 0;
end
try
    suffStat.SXtdfdx;
catch
    suffStat.SXtdfdx = 0;
end
try
    suffStat.SXd2gdx2;
catch
    suffStat.SXd2gdx2 = 0;
end
try
    suffStat.trSx;
catch
    suffStat.trSx = 0;
end
try
    suffStat.SX;
catch
    suffStat.SX = 0;
end
try
    suffStat.SX0;
catch
    suffStat.SX0 = 0;
end

if options.binomial
    try
        suffStat.logL;
    catch
        suffStat.logL = -Inf;
    end
end

%---------------------------------------------------------------------
%------- VBA_Initialize

dim = options.dim;

posterior = options.priors;

posterior.muX = sparse(0,dim.n_t);
indIn = options.params2update.phi;

% watch out about display setting
options.DisplayWin = 0;

%----------------------------------------------------------------------

% Gauss-Newton update of the observation parameters
% !! When the observation function is @VBA_odeLim, this Gauss-Newton update
% actually implements a gradient ascent on the variational energy of the
% equivalent deterministic DCM.

%   for binomial and continuous data

if ~options.OnLine && options.verbose
    fprintf(1,'Deriving prior''s sufficient statistics ...')
    fprintf(1,'%6.2f %%',0)
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
muPhi0 = options.priors.muPhi;
Phi = muPhi0;
Phi(indIn) = posterior.muPhi(indIn);
dphi0 = muPhi0-Phi;
dy = zeros(dim.p,dim.n_t);
vy = zeros(dim.p,dim.n_t);
gx = zeros(dim.p,dim.n_t);


if ~options.binomial
    iQy = options.priors.iQy;
    dy2 = 0;
    Sphid2gdphi2 = 0;
    kernel = zeros(dim.n_phi,dim.n_phi);
    % Get precision parameters
    sigmaHat = posterior.a_sigma./posterior.b_sigma;
    
    
else
    logL = 0;
    
end


if isequal(options.g_fname,@VBA_odeLim)
    muX = zeros(options.inG.old.dim.n,dim.n_t);
    SigmaX = cell(dim.n_t,1);
end
div = 0;

%--- Loop over time series ---%
for t=1:dim.n_t
    
    % evaluate observation function at current mode
    try
        [gx(:,t),dG_dX,dG_dPhi,d2G_dXdPhi] = VBA_evalFun('g',posterior.muX(:,t),Phi,u(:,t),options,dim,t);
        if isweird(gx(:,t))
            disp('')
            disp('Error: could not initialize VB scheme: model generates NaN or Inf!')
            posterior = [];
            return
            
        end
    catch ME
        posterior = [];
        
        name = ME.stack.name; % name of the function in which error occured
        line = ME.stack.line; % number of the line in which error occured
        file = ME.stack.file; % location of the function in which the error occured
        link = ['<a href = "matlab: open ',name,'">',name,'</a>'];
        
        disp('')
        disp(['Error: could not initialize VB scheme : check function "',link,'" at line ', num2str(line)])
        disp('---------------')
        disp('CAUSE :')
        disp([ME.message]); 
        fid = fopen([name,'.m']);
        for l = 1:line
            codeline = fgets(fid);
        end
        fclose(fid);
        disp(codeline)
        disp('---------------')     
        if isequal(ME.message,'Subscripted assignment dimension mismatch.')
            disp('Output dimensions of either the observation or evolution function are incorrect')
        end
        disp(' ')

        return
    end
    
    
    if ~options.binomial
        % mean-field terms
        Sphid2gdphi2 = Sphid2gdphi2 + trace(dG_dPhi*iQy{t}*dG_dPhi'*posterior.SigmaPhi);
        
        % error terms
        dy(:,t) = y(:,t) - gx(:,t);
        dy2 = dy2 + dy(:,t)'*iQy{t}*dy(:,t);
        if dim.n > 0 && ~options.ignoreMF
            A1g = reshape(permute(d2G_dXdPhi,[1,3,2]),dim.p*dim.n,dim.n_phi)';
            A2g = A1g*kron(iQy{t},posterior.SigmaX.current{t});
            kernel = kernel + A2g*A1g';
        end
        
        % Predictive density (data space)
        V = dG_dPhi'*posterior.SigmaPhi*dG_dPhi + (1./sigmaHat).*VB_inv(iQy{t},[]);
        if dim.n > 0
            V = V + dG_dX'*posterior.SigmaX.current{t}*dG_dX;
        end
        vy(:,t) = diag(V);
        
    else
        
        % fix numerical instabilities
        gx(:,t) = checkGX_binomial(gx(:,t));
        
        % predicted variance over binomial data
        vy(:,t) = gx(:,t).*(1-gx(:,t));
        
        % remove irregular trials
        yin = find(~options.isYout(:,t));
        
        % accumulate log-likelihood
        logL = logL + y(yin,t)'*log(gx(yin,t)) + (1-y(yin,t))'*log(1-gx(yin,t));
        
        % prediction error
        dy(yin,t) = y(yin,t) - gx(yin,t);
    end
    
    
    
    % store states dynamics if ODE mode
    if isequal(options.g_fname,@VBA_odeLim)
        % get sufficient statistics of the hidden states from unused i/o in
        % VBA_evalFun.
        muX(:,t) = dG_dX;
        SigmaX{t} = d2G_dXdPhi'*posterior.SigmaPhi*d2G_dXdPhi;
    end
    
    % Display progress
    if mod(t,dim.n_t./10) < 1
        if ~options.OnLine && options.verbose
            fprintf(1,repmat('\b',1,8))
            fprintf(1,'%6.2f %%',100*t/dim.n_t)
        end
    end
    
    % Accelerate divergent update
    if isweird({dy,dG_dPhi,dG_dX})
        div = 1;
        disp(' ')
        disp('Error: could not initialize parameter''s posterior density!')
        disp('Error: Output of evolution or observation function is weird (nan or inf)   ...')
        disp(' ')
        posterior = [];
        break
       
    end
    
end

% Display progress
if ~options.OnLine  && options.verbose % && init
    fprintf(1,repmat('\b',1,8))
    fprintf(' OK.')
    fprintf('\n')
end

% update sufficient statistics
suffStat.Sphi = 0.5*length(indIn)*log(2*pi*exp(1)) + 0.5*VBA_logDet(posterior.SigmaPhi,indIn);
suffStat.gx = gx;
suffStat.dy = dy;
suffStat.vy = vy;
suffStat.dphi = dphi0;
if isequal(options.g_fname,@VBA_odeLim)
    suffStat.muX = muX;
    suffStat.SigmaX = SigmaX;
end
suffStat.div = div;


if ~options.binomial
    suffStat.Sphid2gdphi2 = Sphid2gdphi2;
    suffStat.Sphid2gdphidx = trace(kernel*posterior.SigmaPhi);
    suffStat.dy2 = dy2;
    % Add hyperparameter entropy
    suffStat.Ssigma = gammaln(posterior.a_sigma) ...
        - log(posterior.b_sigma) ...
        + (1-posterior.a_sigma).*psi(posterior.a_sigma) ...
        + posterior.a_sigma;
    
else
    suffStat.logL = logL;
end




end

