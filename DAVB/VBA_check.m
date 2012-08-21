function [options,u,dim] = VBA_check(y,u,f_fname,g_fname,dim,options)
% checks (and fills in default) optional inputs to VB_NLStateSpaceModel.m
% function [options,u,dim] = VBA_check(y,u,f_fname,g_fname,dim,options)
% This function checks the consistency of the arguments to the
% VBA_NLStateSpaceModel.m function. Importantly, it fills in the missing
% priors and optional arguments.
% NB: if the user-specified priors on the stochastic innovations precision
% is a delta Dirac at infinity (priors.a_alpha = Inf and priors.b_alpha =
% 0), the priors are modified to invert an ODE-like state-space model, i.e.
% the equivalent deterministic DCM.


% Fills in dim structure and default input if required
try
    dim.n;
    dim.n_theta;
    dim.n_phi;
catch
    error('Provide dimensions of the model!')
end
try
    dim.n_t;
    dim.p;
catch
    dim.n_t = size(y,2);
    dim.p = size(y,1);
end
if isempty(u)
    u = zeros(1,dim.n_t);
else
    u = full(u);
end
dim.u = size(u,1);



%--------- Check options structure ----------%

% Micro-time resolution
if ~isfield(options,'decim')
    options.decim = 1;
end
% Micro-resolution input
if ~isfield(options,'microU')
    options.microU = 0;
end
% Optional (internal) parameters of the evolution function
if ~isfield(options,'inF')
    options.inF = [];
end
% Optional (internal) parameters of the observation function
if ~isfield(options,'inG')
    options.inG = [];
end
% Check analytical gradients against numerical gradients
if ~isfield(options,'checkGrads')
    options.checkGrads = 0;
end
% Flag for initial conditions update
if ~isfield(options,'updateX0')
    options.updateX0 = 1;
end
% Flag for hyperparameters update
if ~isfield(options,'updateHP')
    options.updateHP = 1;
end
% time lag of the short-sighted backward-pass
if ~isfield(options,'backwardLag')
    options.backwardLag = 1;
end
options.backwardLag = max([floor(options.backwardLag),1]);
% Maximum number of iterations
if ~isfield(options,'MaxIter')
    options.MaxIter = 32;
end
% Maximum number of iterations for the initialization
if ~isfield(options,'MaxIterInit')
    options.MaxIterInit = 32;
end
% Minimum number of iterations
if ~isfield(options,'MinIter')
    options.MinIter = 1;
end
% Minimum change in the free energy
if ~isfield(options,'TolFun')
    options.TolFun = 2e-2;
end
% VB display window
if ~isfield(options,'DisplayWin')
    options.DisplayWin = 1;
end
% Ignore mean-field approximation additional terms (Laplace approx.)?
if ~isfield(options,'ignoreMF')
    options.ignoreMF = 1;
end
% Gauss-Newton ascent on free (1) or variational (0) energy
if ~isfield(options,'gradF')
    options.gradF = 0;
end
% Gauss-Newton coordinate ascent maximum number of iterations:
if ~isfield(options,'GnMaxIter')
    options.GnMaxIter = 32;
end
% Gauss-Newton coordinate ascent threshold on relative variational energy:
if ~isfield(options,'GnTolFun')
    options.GnTolFun = 1e-5;
end
% Gauss-Newton inner-loops figures
if ~isfield(options,'GnFigs')
    options.GnFigs = 0;
end
% matlab window messages
if ~isfield(options,'verbose')
    options.verbose = 1;
end
% On-line version (true when called from VBA_OnLineWrapper.m)
if ~isfield(options,'OnLine')
    options.OnLine = 0;
end
% Free energy calculus at equilibrium?
if ~isfield(options,'Laplace')
    options.Laplace = 1;
end
% delays
if ~isfield(options,'delays')
    options.delays = [];
end
% binomial data
if ~isfield(options,'binomial')
    options.binomial = 0;
end
if options.binomial
    options.ignoreMF = 1;
    if ~isbinary(y(:))
        error('Data should be binary!')
    end
end
if ~isfield(options,'isYout')
    options.isYout = zeros(dim.p,dim.n_t);
end
if ~isfield(options,'skipf')
    options.skipf = zeros(1,dim.n_t);
end
% split-Laplace VB?
if ~isfield(options,'nmog')
    options.nmog = 1;
end
% Name of the display figure
if ~isfield(options,'figName')
    try,aa=options.priors.a_alpha;catch,aa=[];end
    if ~isinf(aa)
        options.figName = 'VB-Laplace identification of stochastic systems';
    else
        if isempty(f_fname)
            options.figName = 'VB-Laplace inversion of static model';
        else
            options.figName = 'VB-Laplace identification of non stochastic systems';
        end
    end
end

% Deal with micro-resolution input
u = VBA_getU(u,options,dim,'2macro');
VBA_disp(' ',options)

%---------- Check priors -----------%
options.params2update.phi = 1:dim.n_phi;
options.params2update.theta = 1:dim.n_theta;
options.params2update.x0 = 1:dim.n;
for t=1:dim.n_t
    options.params2update.x{t} = 1:dim.n;
end
if isfield(options,'priors')
    % Get user-specified priors
    priors = options.priors;
    fn = fieldnames(priors);
    priors0 = VBA_priors(dim,options);
    fn0 = fieldnames(priors0);
    io = ismember(fn0,fn);
    ind = find(io==0);
    if ~isempty(ind)
        VBA_disp('Warning: could not find priors:',options)
        for i = 1:length(ind)
            VBA_disp(['      - ',fn0{ind(i)}],options);
            eval(['priors.',fn0{ind(i)},'=priors0.',fn0{ind(i)},';',])
        end
        VBA_disp('---> Using default(non-informative) priors',options)
    end
    % check dimension and infinite precision priors
    if dim.n_theta > 0 % This finds which evolution params to update
        dpc = diag(priors.SigmaTheta);
        iz = find(dpc==0);
        if ~isempty(iz)
            options.params2update.theta = setdiff(1:dim.n_theta,iz);
        end
    end % This finds which observation params to update
    if dim.n_phi > 0
        dpc = diag(priors.SigmaPhi);
        iz = find(dpc==0);
        if ~isempty(iz)
            options.params2update.phi = setdiff(1:dim.n_phi,iz);
        end
    end  % This finds which initial conditions to update
    if dim.n > 0
        dpc = diag(priors.SigmaX0);
        iz = find(dpc==0);
        if ~isempty(iz)
            options.params2update.x0 = setdiff(1:dim.n,iz);
            if isempty(options.params2update.x0)
                options.updateX0 = 0;
            end
        end
        for t=1:dim.n_t
            dpc = diag(priors.iQx{t});
            iz = find(isinf(dpc));
            if ~isempty(iz)
                options.params2update.x{t} = setdiff(1:dim.n,iz);
            end
        end
    end
else % Build default (non-informative) priors
    priors = VBA_priors(dim,options);
end

% ensure excluded data consistency
if ~options.binomial
    for t=1:dim.n_t
        diQ = diag(priors.iQy{t}).*~options.isYout(:,t);
        options.isYout(:,t) = ~diQ;
        priors.iQy{t} = diag(diQ)*priors.iQy{t}*diag(diQ);
    end
end

% store evolution/observation function handles
options.f_fname = f_fname;
options.g_fname = g_fname;


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

% Special cases check
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
    else
        if ~options.binomial
            % Derive marginalization operators for the lagged Kalman filter
            n = dim.n;
            lag = options.backwardLag + 1;
            options.lagOp.C = [zeros(n,n*(lag-1)),eye(n)];
            options.lagOp.D = [zeros(n,n*(lag-2)),eye(n),zeros(n,n)];
            options.lagOp.E = [eye(n*(lag-1)),zeros(n*(lag-1),n)];
            options.lagOp.Eu = [zeros(n*(lag-1),n),eye(n*(lag-1))];
            options.lagOp.M = [eye(n),zeros(n,n*(lag-1))];
        else
            VBA_disp('Error: stochastic system not supported for binomial data!',options)
        end
    end
end





