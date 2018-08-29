function [f_fname_multi,g_fname_multi,dim_multi,options,u_multi] = VBA_multisession_expand(f_fname,g_fname,dim,options,u)
% [f_fname_multi,g_fname_multi,dim_multi,options,u_multi] = VBA_multisession_expand(f_fname,g_fname,dim,options,u)
% prepare model to be fitted to multisession data. In particular, it:
% - duplicates the parameter set for each session, unless specified 
%   otherwise by multisession.fixed structure
% - wrap the evolution and observation functions so the correct states and
%   parameters are used in each session
% The reverse operation can be performed by VBA_multisession_factor
% to reconstruct within session parameters and states from the expanded model
%% check args

%% check if multisession is required
if ~isfield(options,'multisession') || ~isfield(options.multisession,'split') || numel(options.multisession.split) < 2 ... % no need for multisession
        || (isfield(options.multisession,'expanded') && options.multisession.expanded) % already went there before, no need to expand again
    f_fname_multi = f_fname;
    g_fname_multi = g_fname;
    dim_multi = dim;
    u_multi = u;
    return;
end

    
%% extract sessions
% ensure horizontal vector for looping...
options.multisession.split = options.multisession.split(:)';

if sum(options.multisession.split) ~= dim.n_t
    error('*** Multisession: partition covers %d datapoints but data has %d.',sum(options.multisession.split),dim.n_t);
end

n_session = numel(options.multisession.split);
session_id = ones(1,dim.n_t);
for i=cumsum(options.multisession.split(1:end-1))
    session_id(i+1:end) = session_id(i+1:end) + 1;
end

for i=cumsum(options.multisession.split(1:end-1))
    session_id(i+1) = -session_id(i+1) ; % signal beginning of a session
end

% = append session number to inputs
if options.microU
    session_id = repmat(session_id,options.decim,1);
    session_id = session_id(:)';
end

dim_multi = dim;

u_multi = [u; session_id] ;
dim_multi.u = dim.u+1 ;


%% duplicate parameters
options = VBA_fillInPriors (dim, options);
priors_multi = options.priors;

% = get indexes of duplicated parameters
theta_multi = 1:dim.n_theta;
phi_multi   = 1:dim.n_phi;
X0_multi    = 1:dim.n;

% = restrict fixed parameters
options.multisession = VBA_check_struct(options.multisession,'fixed',struct);
fixed = VBA_check_struct(options.multisession.fixed, ...
    'theta'   , []   , ...
    'phi'     , []   , ...
    'X0'      , []     ...
);
% syntactic sugar handling 
if isequal(fixed.theta ,'all'), fixed.theta = 1:dim.n_theta; end
if isequal(fixed.phi   ,'all'), fixed.phi   = 1:dim.n_phi; end
%if isequal(fixed.X0    ,'all'), fixed.X0    = 1:dim.n; end
    
theta_multi = setdiff(theta_multi,fixed.theta);
phi_multi   = setdiff(phi_multi  ,fixed.phi  );
%X0_multi    = setdiff(X0_multi   ,fixed.X0   );

% = expand (duplicate) priors and dimensions to cover all sessions
priors = options.priors;

try
[priors_multi.muX0, priors_multi.SigmaX0, dim_multi.n] ...
    = expand_param(priors.muX0,priors.SigmaX0,X0_multi,n_session) ;
end
try
[priors_multi.muTheta, priors_multi.SigmaTheta, dim_multi.n_theta] ...
    = expand_param(priors.muTheta,priors.SigmaTheta,theta_multi,n_session) ;
end
[priors_multi.muPhi, priors_multi.SigmaPhi, dim_multi.n_phi] ...
    = expand_param(priors.muPhi,priors.SigmaPhi,phi_multi,n_session) ;


% = restrict initial hidden states
    % enforce covariance across states (duplication is needed for evolution
    % independance)
%     for i = fixed.X0
%         X0_cor = i + (0:n_session-1)*dim.n;
%         priors_multi.SigmaX0(X0_cor,X0_cor) = priors_multi.SigmaX0(i,i);
%     end
    


options.priors = priors_multi;
if isfield(options,'dim')
    options=rmfield(options,'dim');
end

% = precompute param indexes for each session
indices.X0 = param_indices(dim.n,X0_multi,n_session);
indices.theta = param_indices(dim.n_theta,theta_multi,n_session);
indices.phi = param_indices(dim.n_phi,phi_multi,n_session);

multisession.indices = indices;

%% set new evolution and observation functions

multisession.f_fname = f_fname;
multisession.g_fname = g_fname;
multisession.dim = dim_multi;
multisession.X0_multi = X0_multi;
multisession.theta_multi = theta_multi;
multisession.phi_multi = phi_multi;

options.inF = {options.inF multisession} ;
options.inG = {options.inG multisession} ;

f_fname_multi = @f_multi;
g_fname_multi = @g_multi;

options.multisession.expanded = true;
end

%% wrappers for evolution observation functions 

% = wrapper for the evolution function
function  [fx,dF_dX,dF_dTheta] = f_multi(Xt,Theta,ut,in)
    
    % extract options
    inF = in{1};
    multisession = in{2};
    % extract session wise states and params
    session_id = abs(ut(end));
    idx_X0 = multisession.indices.X0(:,session_id);
    idx_theta = multisession.indices.theta(:,session_id);
    
    if ut(end) < 1
        Xt = f_multi(Xt,Theta,[ut(1:end-1);session_id] ,in);
    end
    % call original function
    nout = nargout(multisession.f_fname);
    [output{1:nout}] = multisession.f_fname( ...
    Xt(idx_X0), ...
    Theta(idx_theta), ...
    ut(1:end-1),...
    inF) ;

    % store evolution
    fx = Xt;
    fx(idx_X0) = output{1};
    
      
    % store derivatives if possible
    if nout>=2
        dF_dX = eye(multisession.dim.n,multisession.dim.n);
        dF_dX(idx_X0,idx_X0) = output{2} ;
    else
        dF_dX = [];
    end
    
    if nout>=3
        dF_dTheta = zeros(multisession.dim.n_theta,multisession.dim.n);
        dF_dTheta(idx_theta,idx_X0) = output{3} ;
    else
        dF_dTheta = [];
    end
    
end

% = wrapper for the observation function
function  [gx,dG_dX,dG_Phi] = g_multi(Xt,Phi,ut,in)
    
    % extract options
    inG = in{1};
    multisession = in{2};
    
    % extract session wise states and params
    session_id = abs(ut(end));
    idx_X0 = multisession.indices.X0(:,session_id);
    idx_phi = multisession.indices.phi(:,session_id);
    
    % call original function
    nout = nargout(multisession.g_fname);
    [output{1:nout}] = multisession.g_fname( ...
    Xt(idx_X0), ...
    Phi(idx_phi), ...
    ut(1:end-1),...
    inG) ;

    % store observation
    gx = output{1};
    
    % store derivatives if possible
    if nout>=2
        dG_dX = zeros(multisession.dim.n,numel(gx));
        dG_dX(idx_X0,:) = output{2} ;
    else
        dG_dX = [];
    end
    
    if nout>=3
        dG_Phi = zeros(multisession.dim.n_phi,numel(gx));
        dG_Phi(idx_phi,:) = output{3} ;
    else
        dG_Phi = [];
    end
    
end

%% some shortcuts
function [mu_multi,sigma_multi,dim_multi] = expand_param(mu,sigma,idx,n_session)
% concatenate priors means and variances accross sessions

% = initial and ifnal dimensions
n1 = numel(mu);
n2 = n1 + (n_session-1)*numel(idx);

% = means
mu = mu(:);
mu_multi = [mu ; repmat(mu(idx),n_session-1,1)];

% = variances
sigma_temp = kron(eye(n_session-1), sigma(idx,idx));
sigma_multi = [sigma,           zeros(n1,n2-n1);
               zeros(n2-n1,n1)  sigma_temp];

% = dimension          
dim_multi = n2;                    
                    
end

function indices = param_indices(n,idx,n_session)
% return for each session (column) the indices of n parameters in the 
% expanded priorsgiven only those in idx are duplicated

indices = repmat(1:n,n_session,1)';
for k=2:n_session
    indices(idx,k) = n + (k-2)*numel(idx) + (1:numel(idx));
end

end

