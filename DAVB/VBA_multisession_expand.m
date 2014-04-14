function [f_fname_multi,g_fname_multi,dim_multi,options,u_multi] = VBA_multisession_expand(f_fname,g_fname,dim,options,u)


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
if sum(options.multisession.split) ~= dim.n_t
    error('*** Multisession: partition covers %d datapoints but data has %d.',sum(options.multisession.split),dim.n_t);
end

n_session = numel(options.multisession.split);
session_id = ones(1,dim.n_t);
for i=cumsum(options.multisession.split(1:end-1))
    session_id(i+1:end) = session_id(i+1:end) + 1;
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
priors_multi = options.priors;

if ~isfield(options.multisession,'fixed')
    options.multisession.fixed = struct();
end

if ~isfield(options.multisession.fixed,'X0') 
    options.multisession.fixed.X0 = [];
end
if ~isfield(options.multisession.fixed,'theta') 
    options.multisession.fixed.theta = [];
end
if ~isfield(options.multisession.fixed,'phi') 
    options.multisession.fixed.phi = [];
end

% = expand (duplicate) priors and dimensions to cover all sessions
priors = options.priors;

if dim.n > 0
[priors_multi.muX0, priors_multi.SigmaX0, dim_multi.n] = expand_param(priors.muX0,priors.SigmaX0,n_session,options.multisession.fixed.X0) ; 
end
if dim.n_theta > 0
[priors_multi.muTheta, priors_multi.SigmaTheta, dim_multi.n_theta] = expand_param(priors.muTheta,priors.SigmaTheta,n_session,options.multisession.fixed.theta) ;    
end
if dim.n_phi > 0
[priors_multi.muPhi, priors_multi.SigmaPhi, dim_multi.n_phi] = expand_param(priors.muPhi,priors.SigmaPhi,n_session,options.multisession.fixed.phi) ;
end


options.priors = priors_multi;
if isfield(options,'dim')
    options=rmfield(options,'dim');
end

% = precompute param indexes for each session
indices.X0 = param_indices(dim.n,n_session);
indices.theta = param_indices(dim.n_theta,n_session);
indices.phi = param_indices(dim.n_phi,n_session);

multisession.indices = indices;

%% set new evolution and observation functions

multisession.f_fname = f_fname;
multisession.g_fname = g_fname;
multisession.dim = dim_multi;

options.inF = {options.inF multisession} ;
options.inG = {options.inG multisession} ;

f_fname_multi = @f_multi;
g_fname_multi = @g_multi;

options.multisession.expanded = true;
options.multisession.indices = indices;
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
    [output{1:nout}] = feval(multisession.f_fname, ...
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
    [output{1:nout}] = feval(multisession.g_fname, ...
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
function [mu_multi,sigma_multi,dim_multi] = expand_param(mu,sigma,n_session,fixed)
% concatenate priors means and variances accross sessions

% = dimension   
n = numel(mu);
dim_multi =  n*n_session;  

% = means
mu_multi = repmat(mu(:),n_session,1);

% = variances
sigma_multi = kron(eye(n_session),sigma);

for i=fixed
        icor = i + (0:n_session-1)*n;
        sigma_multi(icor,icor) = sigma(i,i);
end          
                    
end

function indices = param_indices(n,n_session)
% return for each session (column) the indices of n parameters in the 
% expanded priorsgiven only those in idx are duplicated

indices = reshape(1:(n*n_session),n,n_session);


end

