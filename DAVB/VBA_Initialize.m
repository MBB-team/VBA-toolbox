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
    % free energy at the prior pdf while testing if their is any errors in
    % model structure, inputs, or outputs
%{    
%         
%     posterior.muX = sparse(0,dim.n_t);
%     indIn = options.params2update.phi;
%     opt = options;
%     suffStat = VBA_getSuffStat(opt);    
%     
%     % watch out about display setting
%     opt.DisplayWin = 0;
%              
%     %----------------------------------------------------------------------
%    
%     % Gauss-Newton update of the observation parameters
%     % !! When the observation function is @VBA_odeLim, this Gauss-Newton update
%     % actually implements a gradient ascent on the variational energy of the
%     % equivalent deterministic DCM.
%     
%     %   for binomial and continuous data 
%     
%     % check if called during initialization
%     if isequal(suffStat,VBA_getSuffStat(opt))
%         init = 1;
%         if ~opt.OnLine && opt.verbose
%             fprintf(1,'Deriving prior''s sufficient statistics ...')
%             fprintf(1,'%6.2f %%',0)
%         end
%     else
%         init = 0;
%     end
%     
%     if opt.DisplayWin && ~init % Display progress
%         if isequal(opt.g_fname,@VBA_odeLim)
%             STR = 'VB Gauss-Newton on observation/evolution parameters... ';
%         else
%             STR = 'VB Gauss-Newton on observation parameters... ';
%         end
%         set(opt.display.hm(1),'string',STR);
%         set(opt.display.hm(2),'string','0%');
%         drawnow
%     end
%     
%     % Clear persistent variables if ODE mode
%     if isequal(opt.g_fname,@VBA_odeLim) || ...
%             isequal(opt.g_fname,@VBA_smoothNLSS)
%         clear VBA_odeLim
%         clear VBA_smoothNLSS
%     end
%     
%     %  Look-up which evolution parameter to update
%     indIn = opt.params2update.phi;
%     
%     % Preallocate intermediate variables
%     muPhi0 = opt.priors.muPhi;
%     Phi = muPhi0;
%     Phi(indIn) = posterior.muPhi(indIn);
%     dphi0 = muPhi0-Phi;
%     dy = zeros(dim.p,dim.n_t);
%     vy = zeros(dim.p,dim.n_t);
%     gx = zeros(dim.p,dim.n_t);
%     
%     
%     if ~opt.binomial
%         iQy = opt.priors.iQy;
%         dy2 = 0;
%         Sphid2gdphi2 = 0;
%         kernel = zeros(dim.n_phi,dim.n_phi);
%         % Get precision parameters
%         sigmaHat = posterior.a_sigma./posterior.b_sigma;
%         
%         
%     else
%         logL = 0;
%         
%     end
%     
%     
%     if isequal(opt.g_fname,@VBA_odeLim)
%         muX = zeros(opt.inG.old.dim.n,dim.n_t);
%         SigmaX = cell(dim.n_t,1);
%     end
%     div = 0;
%     
%     %--- Loop over time series ---%
%     for t=1:dim.n_t
%         
%         % evaluate observation function at current mode
%         try
%             [gx(:,t),dG_dX,dG_dPhi,d2G_dXdPhi] = VBA_evalFun('g',posterior.muX(:,t),Phi,u(:,t),opt,dim,t);
%             if isweird(gx(:,t))
%                 disp('')
%                 disp('Error: could not initialize VB scheme: model generates NaN or Inf!')
%                 posterior = [];
% 
%             end
%         catch ME
%             disp('')
%             disp('Error: could not initialize VB scheme: check your model functions!')
%             disp('- check model dimensions (hidden states and parameters)');
%             disp('- check indices in the function');
%             
%             disp('----------------')
%             disp([ME.getReport]);
%             disp('----------------')
%             
% 
%             posterior = [];
%             return
%         end
% 
%         
%         if ~opt.binomial
%             % mean-field terms
%             Sphid2gdphi2 = Sphid2gdphi2 + trace(dG_dPhi*iQy{t}*dG_dPhi'*posterior.SigmaPhi);
%             
%             % error terms
%             dy(:,t) = y(:,t) - gx(:,t);
%             dy2 = dy2 + dy(:,t)'*iQy{t}*dy(:,t);
%             if dim.n > 0 && ~opt.ignoreMF
%                 A1g = reshape(permute(d2G_dXdPhi,[1,3,2]),dim.p*dim.n,dim.n_phi)';
%                 A2g = A1g*kron(iQy{t},posterior.SigmaX.current{t});
%                 kernel = kernel + A2g*A1g';
%             end
%             
%             % Predictive density (data space)
%             V = dG_dPhi'*posterior.SigmaPhi*dG_dPhi + (1./sigmaHat).*VB_inv(iQy{t},[]);
%             if dim.n > 0
%                 V = V + dG_dX'*posterior.SigmaX.current{t}*dG_dX;
%             end
%             vy(:,t) = diag(V);
%             
%         else
%             
%             % fix numerical instabilities
%             gx(:,t) = checkGX_binomial(gx(:,t));
%             
%             % predicted variance over binomial data
%             vy(:,t) = gx(:,t).*(1-gx(:,t));
%             
%             % remove irregular trials
%             yin = find(~opt.isYout(:,t));
%             
%             % accumulate log-likelihood
%             logL = logL + y(yin,t)'*log(gx(yin,t)) + (1-y(yin,t))'*log(1-gx(yin,t));
%             
%             % prediction error
%             dy(yin,t) = y(yin,t) - gx(yin,t);
%         end
%         
%         
%         
%         % store states dynamics if ODE mode
%         if isequal(opt.g_fname,@VBA_odeLim)
%             % get sufficient statistics of the hidden states from unused i/o in
%             % VBA_evalFun.
%             muX(:,t) = dG_dX;
%             SigmaX{t} = d2G_dXdPhi'*posterior.SigmaPhi*d2G_dXdPhi;
%         end
%         
%         % Display progress
%         if mod(t,dim.n_t./10) < 1
%             if ~init && opt.DisplayWin
%                 set(opt.display.hm(2),'string',[num2str(floor(100*t/dim.n_t)),'%']);
%                 drawnow
%             end
%             if init && ~opt.OnLine && opt.verbose
%                 fprintf(1,repmat('\b',1,8))
%                 fprintf(1,'%6.2f %%',100*t/dim.n_t)
%             end
%         end
%         
%         % Accelerate divergent update
%         if isweird({dy,dG_dPhi,dG_dX})
%             div = 1;
%             break
%         end
%         
%     end
%     
%     % Display progress
%     if ~init && opt.DisplayWin
%         set(opt.display.hm(2),'string','OK');
%         drawnow
%     end
%     if init &&  ~opt.OnLine  && opt.verbose
%         fprintf(1,repmat('\b',1,8))
%         fprintf(' OK.')
%         fprintf('\n')
%     end
%     
%     % update sufficient statistics
%     suffStat.Sphi = 0.5*length(indIn)*log(2*pi*exp(1)) + 0.5*VBA_logDet(posterior.SigmaPhi,indIn);
%     suffStat.gx = gx;
%     suffStat.dy = dy;
%     suffStat.vy = vy;
%     suffStat.dphi = dphi0;
%     if isequal(opt.g_fname,@VBA_odeLim)
%         suffStat.muX = muX;
%         suffStat.SigmaX = SigmaX;
%     end
%     suffStat.div = div;
%     
%     
%     if ~opt.binomial
%         suffStat.Sphid2gdphi2 = Sphid2gdphi2;
%         suffStat.Sphid2gdphidx = trace(kernel*posterior.SigmaPhi);
%         suffStat.dy2 = dy2;
%         % Add hyperparameter entropy
%         suffStat.Ssigma = gammaln(posterior.a_sigma) ...
%             - log(posterior.b_sigma) ...
%             + (1-posterior.a_sigma).*psi(posterior.a_sigma) ...
%             + posterior.a_sigma;
%         
%     else
%         suffStat.logL = logL;
%     end
    
    %----------------------------------------------------------------------
%}    
    [ suffStat, posterior ] = VBA_check_errors(y,u, options);    
    options.init.posterior = posterior;
    options.init.suffStat = suffStat;
    
end

options.init.F = VBA_FreeEnergy(posterior,suffStat,options);



