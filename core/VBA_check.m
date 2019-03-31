function [options,u,dim] = VBA_check(y,u,f_fname,g_fname,dim,options)
% VBA_CHECK 
% [options,u,dim] = VBA_check(y,u,f_fname,g_fname,dim,options)
% Checks (and fills in default) optional inputs to VBA_NLStateSpaceModel.m
% function [options,u,dim] = VBA_check(y,u,f_fname,g_fname,dim,options)
% This function checks the consistency of the arguments to the
% VBA_NLStateSpaceModel.m function. Importantly, it fills in the missing
% priors and optional arguments.
% NB: if the user-specified priors on the stochastic innovations precision
% is a delta Dirac at infinity (priors.a_alpha = Inf and priors.b_alpha =
% 0), the priors are modified to invert an ODE-like state-space model, i.e.
% the equivalent deterministic system.

%% ________________________________________________________________________
%  check model dimension

% set default if static model
if isempty(f_fname)
    dim = VBA_check_struct(dim, ...
        'n_theta', 0, ...
        'n', 0);
end
    
% check if hidden state dimension:
assert( ...
     isfield(dim,'n')          ... exists,
  && isscalar(dim.n)           ... is a scalar value
  && mod(dim.n,1) == 0         ... that is integer
  && dim.n >= 0                ... and non negative
, '*** Please provide the dimension of the hidden state in dim.n');

% check if evolution parameters dimension:
assert( ...
     isfield(dim,'n_theta')    ... exists
  && isscalar(dim.n_theta)     ... is a scalar value
  && mod(dim.n_theta,1) == 0   ... that is integer
  && dim.n_theta >= 0          ... and non negative
, '*** Please provide the dimension of the evolution parameters in dim.n_theta');

% check if observation parameters dimension:
assert( ...
     isfield(dim,'n_phi')      ... exists
  && isscalar(dim.n_phi)       ... is a scalar value
  && mod(dim.n_phi,1) == 0     ... that is integer
  && dim.n_phi >= 0            ... and non negative
, '*** Please provide the dimension of the observation parameters in dim.n_phi');


dim = VBA_check_struct(dim, ...
    'n_t'   , size(y,2), ... % number of trials or time samples
    'p'     , size(y,1)  ... % data dimension
);

if isempty(u)
    u = zeros(0,dim.n_t);
else
    u = full(u);
end
dim.u = size(u,1);

%% ________________________________________________________________________
%  check VBA's options structure

if isfield(options,'binomial') % check for legacy binomial field
    options.sources = struct('type', options.binomial , ...
                             'out' , 1:dim.p          );
    options = rmfield(options, 'binomial');                    
    fprintf('*** WARNING: options.binomial is now deprecated. Please use ''options.sources.type = 1'' to specify a binary observation distribution.\n');
end

% set defaults 
options = VBA_check_struct(options, ...
    'decim'      , 1     , ...     % Micro-time resolution
    'microU'     , 0     , ...     % Micro-resolution input
    'inF'        , []    , ...     % Optional (internal) parameters of the evolution function
    'inG'        , []    , ...     % Optional (internal) parameters of the observation function
    'checkGrads' , 0     , ...     % Check analytical gradients against numerical gradients
    'updateX0'   , 1     , ...     % Flag for initial conditions update
    'updateHP'   , 1     , ...     % Flag for hyperparameters update
    'initHP'     , 1     , ...     % Flag for hyperparameters initialization
    'backwardLag', 1     , ...     % time lag of the short-sighted backward-pass
    'MaxIter'    , 32    , ...     % Maximum number of iterations
    'MaxIterInit', 32    , ...     % Maximum number of iterations for the initialization
    'MinIter'    , 0     , ...     % Minimum number of iterations
    'TolFun'     , 2e-2  , ...     % Minimum change in the free energy
    'DisplayWin' , 1     , ...     % VB display window
    'gradF'      , 0     , ...     % Gauss-Newton ascent on free (1) or variational (0) energy
    'GnMaxIter'  , 32    , ...     % Gauss-Newton coordinate ascent maximum number of iterations:
    'GnTolFun'   , 1e-5  , ...     % Gauss-Newton coordinate ascent threshold on relative variational energy:
    'GnFigs'     , 0     , ...     % Gauss-Newton inner-loops figures
    'verbose'    , 1     , ...     % matlab window messages
    'OnLine'     , 0     , ...     % On-line version (true when called from VBA_OnLineWrapper.m)
    'delays'     , []    , ...     % delays
    'kernelSize' , 16    , ...     % max lag of Volterra kernels
    'nmog'       , 1     , ...     % split-Laplace VB?
    'UNL'        , 0     , ...     % un-normalized likelihood?
    'UNL_width'  , 4     , ...     % for partition function estimation
    'UNL_ng'     , 64      ...     % for partition function estimation
) ;
 
options = VBA_check_struct(options, ...
    'isYout'    , zeros(dim.p,dim.n_t)    , ... % excluded data
    'skipf'     , zeros(1,dim.n_t)        , ... % steady evolution
    'sources'   , struct()                  ...
) ;

options.sources = VBA_check_struct(options.sources, ...
    'type', 0        , ... % default to gaussian
    'out' , 1:dim.p  );
                     

% % retrocompatibility
% if ~isfield(options,'binomial') && numel(options.sources)==1
%         options.binomial = options.sources.type;
% end

    
    
options.backwardLag = min([max([floor(round(options.backwardLag)),1]),dim.n_t]);
options.kernelSize  = min([dim.n_t,options.kernelSize]);

for i=1:numel(options.sources)
    if options.sources(i).type ~= 0 % if binomial
        if ~isempty(y)
            isYoutSource = options.isYout(options.sources(i).out,:);
            ySource = y(options.sources(i).out,:);
            if ~ VBA_isBinary (ySource(~ isYoutSource))
                error('*** Data should be binary!')
            end
            if islogical(ySource(~isYoutSource))
                error('*** Binary data should be numerical (0,1) not logical (true,false).\n*** Use y=+y or y=double(y) to convert your observations.',[])
            end
        end
    end
end

% split in multisession
[f_fname,g_fname,dim,options,u] = VBA_multisession_expand(f_fname,g_fname,dim,options,u);

% Deal with micro-resolution input
u = VBA_getU(u,options,dim,'2macro');


%% ________________________________________________________________________
%  check priors
[options,options.params2update] = VBA_fillInPriors(dim,options);
priors = options.priors;
% check consistency of priors and model dimensions
assert(size(priors.muX0,1)==dim.n,'*** Dimension of options.priors.muX0 does not match dim.n!')
assert(isequal(size(priors.SigmaX0),dim.n*[1,1]),'*** Dimension of options.priors.SigmaX0 does not match dim.n!')
assert(size(priors.muPhi,1)==dim.n_phi,'*** Dimension of options.priors.muPhi does not match dim.n_phi!')
assert(isequal(size(priors.SigmaPhi),dim.n_phi*[1,1]),'*** Dimension of options.priors.SigmaPhi does not match dim.n_phi!')
assert(size(priors.muTheta,1)==dim.n_theta,'*** Dimension of options.priors.muTheta does not match dim.n_theta!')
assert(isequal(size(priors.SigmaTheta),dim.n_theta*[1,1]),'*** Dimension of options.priors.SigmaTheta does not match dim.n_theta!')

if dim.n>0 && isempty(options.params2update.x0)
    options.updateX0 = 0;
end

% Name of the display figure
if ~isfield(options,'figName')
    if dim.n > 0
        if isinf(priors.a_alpha) && isequal(priors.b_alpha,0)
            options.figName = 'VB-Laplace identification of deterministic system';
        else
            options.figName = 'VB-Laplace identification of stochastic system';
        end
    else
        options.figName = 'VB-Laplace inversion of static model';
    end
end

% % ensure excluded data consistency (gaussian sources)
% TODO: remove once VBA_Iphi deal with isYout
if numel(options.sources) == 1
gsi = find([options.sources.type]==0);
for i=1:numel(gsi)
    for t=1:dim.n_t
        diQ = (diag(priors.iQy{t,i}) == 0) | options.isYout(options.sources(gsi(i)).out,t);
        options.isYout(options.sources(gsi(i)).out,t) = +diQ;
        priors.iQy{t,i} = diag(~diQ)*priors.iQy{t,i}*diag(~diQ);
    end
end
end


% store evolution/observation function handles
options.f_fname = f_fname;
options.g_fname = g_fname;


%% ________________________________________________________________________
%  check inversion specific cases

% split-MoG sufficient statistics
if options.nmog > 1
    np = length(options.params2update.phi);
    [m0,s0,w0] = getMoG4N01(options.nmog,1,0);
    tic,[C] = get_nkdraws(options.nmog,np,0);toc
    nd = size(C,2);
    split.w = zeros(nd,1);
    split.m = zeros(np,nd);
    split.s = zeros(np,nd);
    for i=1:nd
        split.w(i) = prod(w0(C(:,i)));
        split.m(:,i) = m0(C(:,i));
        split.s(:,i) = s0(C(:,i));
    end
    split.w = split.w./sum(split.w);
    options.split = split;
end

% Delay embedding
if dim.n > 0 && sum(options.delays(:)) > 0
    % modify options:
    dim.n_embed = dim.n*max(options.delays(:));
    inF.dim = dim;
    inF.options = options;
    inF.f_fname = f_fname;
    inF.f_nout = nargout(f_fname);
    inG.dim = dim;
    inG.options = options;
    inG.g_fname = g_fname;
    inG.g_nout = nargout(g_fname);
    options.inF = inF;
    options.inG = inG;
    options.f_fname = @f_embed;
    options.g_fname = @g_embed;
    options.delays = [];
    % modify priors:
    Zeros = zeros(inF.dim.n,dim.n_embed);
    for t=1:dim.n_t
        tmp = 1e2*kron(eye(max(inG.options.delays(:))),priors.iQx{t});
        priors.iQx{t} = [   priors.iQx{t}	Zeros
                            Zeros'          tmp    ];
        priors.iQx{t}(isnan(priors.iQx{t})) = 0;
        dpc = diag(priors.iQx{t});
        iz = find(isinf(dpc));
        if ~isempty(iz)
            options.params2update.x{t} = setdiff(1:dim.n,iz);
        else
            options.params2update.x{t} = 1:dim.n;
        end
    end
    priors.muX0 = repmat(priors.muX0,max(inF.options.delays(:))+1,1);
    priors.SigmaX0 = kron(eye(max(inF.options.delays(:))+1),priors.SigmaX0);
    dpc = diag(priors.SigmaX0);
    iz = find(dpc==0);
    if ~isempty(iz)
        options.params2update.x0 = setdiff(1:length(priors.SigmaX0),iz);
    else
        options.params2update.x0 = 1:length(priors.SigmaX0);
    end
    dim.n = dim.n*(max(inF.options.delays(:))+1);
end

% Add other inputs in the options structure:
options.g_nout = nargout(options.g_fname);
options.priors = priors;
options.dim = dim;

% Special case: DCM check
if isequal(options.f_fname,@f_DCMwHRF) && isequal(options.g_fname,@g_HRF3)
    % DCM for fMRI
    [options] = VBA_check4DCM(options);
end

if dim.n > 0
    options.f_nout = nargout(options.f_fname);
    if isinf(priors.a_alpha) && isequal(priors.b_alpha,0)
        % Check whether SDE --> ODE [ p(alpha) = Dirac at 0 ]
        % This will allow the algorithm to shortcut the VB-Laplace extended Kalman
        % smoother, treating the model as non-stochastic DCM. In this version, the
        % hidden states are assumed to obey an ODE, which is the deterministic
        % variant of the more general SDE.
        priors0 = priors;
        dim0 = dim;
        options0 = options;
        % Embed evolution parameters and initial conditions in observation
        % parameters  --> modify prior structure accordingly
        priors.muTheta = [];
        priors.SigmaTheta = [];
        priors.muX0 = [];
        priors.SigmaX0 = [];
        if dim0.n_theta > 0
            priors.muPhi = [priors0.muPhi;priors0.muTheta];
            Zeros = zeros(dim0.n_phi,dim0.n_theta);
            priors.SigmaPhi = [ priors0.SigmaPhi  Zeros
                                Zeros'            priors0.SigmaTheta ];
            options.params2update.phi = [options0.params2update.phi,options0.params2update.theta + dim0.n_phi];
        end
        if options.updateX0
            priors.muPhi = [priors.muPhi;priors0.muX0];
            Zeros = zeros(dim0.n_phi+dim0.n_theta,dim0.n);
            priors.SigmaPhi = [ priors.SigmaPhi  Zeros
                                Zeros'           priors0.SigmaX0 ];
            options.params2update.phi = [options.params2update.phi,options0.params2update.x0 + dim0.n_phi + dim0.n_theta];
        end
        % Modify dim and options structures
        dim.n_theta = 0;
        dim.n_phi = size(priors.muPhi,1);
        dim.n = 0;
        options.f_fname = [];
        options.f_nout = 0;
        options.g_fname = @VBA_odeLim;
        options.g_nout = nargout(@VBA_odeLim);
        options.inF = [];
        options.inG = [];
        options.inG.old.dim = dim0;
        options.inG.old.options = options0;
        options.inG.dim = dim;
        options.priors = priors;
        options.dim = dim;
    else % stochastic inversion
        % Derive marginalization operators for the lagged Kalman filter
        n = dim.n;
        lag = options.backwardLag + 1;
        options.lagOp.C = [zeros(n,n*(lag-1)),eye(n)];
        options.lagOp.D = [zeros(n,n*(lag-2)),eye(n),zeros(n,n)];
        options.lagOp.E = [eye(n*(lag-1)),zeros(n*(lag-1),n)];
        options.lagOp.Eu = [zeros(n*(lag-1),n),eye(n*(lag-1))];
        options.lagOp.M = [eye(n),zeros(n,n*(lag-1))];
    end
end





