function [posterior,suffStat,options] = VBA_Initialize(y,u,f_fname,g_fname,dim,options)
% initializes approximate posterior of sDCM unknown variables
% function [posterior,suffStat] = VBA_Initialize(y,u,muX0,f_fname,g_fname,options)
%
% This function initializes the parameters and hidden-states which have to
% be estimated using the VBA approach to NL state-space model.
% NB: stochastic DCMs are initialized with their deterministic limit.
% This ensures that stochastic innovations (state noise) only explain what
% could not be explained by the deterministic DCM.

% Initialize posterior with priors
posterior = options.priors;
posterior = rmfield(posterior,'iQy');
try, posterior = rmfield(posterior,{'iQx','AR'}); end

if dim.n > 0        % if any hidden-states to be estimated
    
    % Initialize with VB-Laplace for non-stochastic systems:
    % Get approximate posterior for evolution/observation parameters from
    % ODE limit to the state-space model.
    if dim.n_phi + dim.n_theta > 0
        
        options0 = options;
        options0.figName = 'VB-Laplace initialization: non stochastic system';
        options0.priors.a_alpha = Inf;
        options0.priors.b_alpha = 0;
        options0.MaxIter = options.MaxIterInit;

        if ~options.priors.AR
            u = VBA_getU(u,options0,dim,'back2micro');
            [in.out.options,u,in.out.dim] = VBA_check(...
                y,...
                u,...
                options0.f_fname,...
                options0.g_fname,...
                dim,...
                options0);
        else % AR(1) state noise
            % VBA_check has already done most of the job...
            options0.f_fname = [];
            options0.f_nout = [];
            options0.inF = [];
            options0.g_fname = @VBA_odeLim;
            options0.g_nout = nargout(@VBA_odeLim);
            in.out.options = options0;
            in.out.dim = dim;
            in.out.dim.n = 0;
        end
        
        [in.posterior,in.out.suffStat,in.out.options] = VBA_Initialize(...
            y,...
            u,...
            f_fname,...
            g_fname,...
            in.out.dim,...
            in.out.options);
        
        in.out.F = VBA_FreeEnergy(in.posterior,in.out.suffStat,in.out.options);
        in.out.u = u;
        
        [posterior,out] = VBA_NLStateSpaceModel(y,[],[],[],[],[],in);
        
        % Store deterministic inversion
        options.init.posterior = posterior;
        options.init.out = out;
        
        if options.init0
            posterior = options.priors;
            % Get hidden states prior predictive density
            [posterior.muX,posterior.SigmaX.current,suffStat] = VBA_EKF(y,u,posterior,dim,options,2);
            % NB: deterministic fit can still be compared to the stochastic
            % variant. It is just not used as an initialization...
        else
            % Only reinitialize precision hyperparameters
            posterior.a_alpha = options.priors.a_alpha;
            posterior.b_alpha = options.priors.b_alpha;
            if ~options.binomial
                posterior.a_sigma = options.priors.a_sigma;
                posterior.b_sigma = options.priors.b_sigma;
            end
            suffStat = out.suffStat;
        end

        In = eye(dim.n);
        alphaHat = options.priors.a_alpha./options.priors.b_alpha;
        if ~~options.priors.AR
            posterior.muX = zeros(dim.n,dim.n_t);
            for t = 1:dim.n_t
                posterior.SigmaX.current{t} = 1./alphaHat*In;
            end
        else
            for t = 1:dim.n_t
                posterior.SigmaX.current{t} = posterior.SigmaX.current{t} + 1./alphaHat*In;
            end
        end
        
    else    % if no evolution/observation parameters to be estimated
        
        % Get hidden states prior predictive density
        [posterior.muX,posterior.SigmaX.current,suffStat] = VBA_EKF(y,u,posterior,dim,options,2);
        
        % Emulate deterministic inversion
        options.init.posterior = posterior;
        options.priors.muX = posterior.muX;
        options.priors.SigmaX = posterior.SigmaX;
        options.init.out = struct(...
            'options',options,...
            'u',VBA_getU(u,options,dim,'back2micro'),...
            'y',y,...
            'dim',dim,...
            'it',0,...
            'suffStat',suffStat,...
            'date',clock,...
            'dt',0);
        
        [suffStat] = VBA_getSuffStat(options,suffStat);
        options.init.out.F = VBA_FreeEnergy(posterior,suffStat,options);
        
    end
    
    posterior.SigmaX.inter = cell(dim.n_t-1,1);
    In = speye(dim.n);
    for t = 1:dim.n_t-1
        posterior.SigmaX.inter{t} = 0*full(In);
    end
    options.priors.muX = posterior.muX;
    options.priors.SigmaX = posterior.SigmaX;

    % Add hyperparameter entropies
    if ~options.binomial
        suffStat.Ssigma = gammaln(posterior.a_sigma) ...
            - log(posterior.b_sigma) ...
            + (1-posterior.a_sigma).*psi(posterior.a_sigma) ...
            + posterior.a_sigma;
    end
    suffStat.Salpha = gammaln(posterior.a_alpha) ...
        - log(posterior.b_alpha) ...
        + (1-posterior.a_alpha).*psi(posterior.a_alpha) ...
        + posterior.a_alpha;
    
    % Get sufficient statistics and free energy
    opt0 = options;
    opt0.DisplayWin = 0;
    if ~options.binomial
        [o1,o2,o3,suffStat] = VBA_IX_lagged(...
            posterior.muX,...
            y,...
            posterior,...
            suffStat,...
            dim,...
            u,...
            opt0);
    else
        [o1,o2,o3,suffStat] = VBA_IX_binomial(...
            posterior.muX,...
            y,...
            posterior,...
            suffStat,...
            dim,...
            u,...
            opt0);
    end
    
    
else        % if no hidden states, initialize observation parameters
    

    % This fills in the sufficient statistics structure, to evaluate the
    % free energy at the prior pdf.
    posterior.muX = sparse(0,dim.n_t);
    indIn = options.params2update.phi;
    opt = options;
    suffStat = VBA_getSuffStat(opt);
    opt.DisplayWin = 0;
    if ~options.binomial
        % Get sufficient statistics for the observation parameters
        [o1,o2,o3,suffStat] = VBA_Iphi(...
            posterior.muPhi(indIn),...
            y,...
            posterior,...
            suffStat,...
            dim,...
            u,...
            opt);
        % Add hyperparameter entropy
        suffStat.Ssigma = gammaln(posterior.a_sigma) ...
            - log(posterior.b_sigma) ...
            + (1-posterior.a_sigma).*psi(posterior.a_sigma) ...
            + posterior.a_sigma;
    else
        [o1,o2,o3,suffStat] = VBA_Iphi_binomial(...
            posterior.muPhi(indIn),...
            y,...
            posterior,...
            suffStat,...
            dim,...
            u,...
            opt);
    end
    
    options.init.posterior = posterior;
    options.init.suffStat = suffStat;
    
end

if suffStat.div
    disp(' ')
    disp('Error: could not initialize parameter''s posterior density!')
    disp(' ')
else
    % Fill in potential missing sufficient statistics
    [suffStat] = VBA_getSuffStat(options,suffStat);
    options.init.F = VBA_FreeEnergy(posterior,suffStat,options);
end





