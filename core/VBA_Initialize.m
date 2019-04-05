function [posterior,suffStat,options] = VBA_Initialize(y,u,dim,options)
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

if dim.n > 0
    
    % Initialize with VB-Laplace for non-stochastic systems:
    % Get approximate posterior for evolution/observation parameters from
    % ODE limit to the state-space model.
    
    options0 = options;
    options0.figName = 'VB-Laplace initialization: deterministic system';
    options0.priors.a_alpha = Inf;
    options0.priors.b_alpha = 0;
    options0.MaxIter = options.MaxIterInit;
    u = VBA_getU(u,options0,dim,'back2micro');
    [posterior,out] = VBA_NLStateSpaceModel(y,u,options.f_fname,options.g_fname,dim,options0);
    
    % Store deterministic inversion
    options.init.posterior = posterior;
    options.init.out = out;
    options.priors.muX = out.options.priors.muX;
    options.priors.SigmaX = out.options.priors.SigmaX;
    suffStat = out.suffStat;
    
    % Reinitialize precision hyperparameters posterior and sufficient stats
    posterior.a_alpha = options.priors.a_alpha;
    posterior.b_alpha = options.priors.b_alpha;

   if any([options.sources.type]==0) && ~options.initHP
       posterior.a_sigma = options.priors.a_sigma;
       posterior.b_sigma = options.priors.b_sigma;
   end
   
   % try to pass on display handles to main inversion
   try
       options.display = out.options.display;
   end

else
    
    [ suffStat, posterior ] = VBA_check_errors(y,u,options);
    options.init.posterior = posterior;
    options.init.suffStat = suffStat;
    
end

options.init.F = VBA_FreeEnergy(posterior,suffStat,options);






