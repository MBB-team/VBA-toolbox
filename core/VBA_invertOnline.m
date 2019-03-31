function [posterior, out] = VBA_invertOnline(y, u, f_fname, g_fname, dim, options)
% This function inverts a nonlinear state-space model 'on line', i.e.
% estimate hidden states, observation/evolution parameters and
% hyperparameters in 'real-time', by adding the time samples one after the
% other. It returns a time-dependent posterior density over all unknown
% variables. See VBA_NLStateSpaceModel.m for I/O arguments.
% NB: this online wrapper does not deal with ODE state-space models.


for i = 1:size(y,2)
    options.priors.iQy{i,1} = 2; 
end

options.OnLine = 1;

if ~ isfield(options, 'isYout')
    isYout = zeros(size(y));
else
    isYout = options.isYout;
end
%
options_online = options;

figure;
h1 = subplot(2,1,1);
h2 = subplot(2,1,2);

in = NaN;
for t = 1 : size(y,2)
    
    t
    options_online.updateX0 = +(t == 1);
    
    options_online.isYout = ones(size(y));
    options_online.isYout(:,t) = isYout(:,t);
    
    
    [posterior(t),out(t)] = VBA_NLStateSpaceModel(y,u,f_fname,g_fname,dim,options_online, in);

    %[posterior,suffStat,options] = updatePSS(...
    %       OL_posterior,OL_out,posterior,suffStat,1,options);
    
    
    plot(h1,1:t,[posterior.muPhi]);
    plot(h2,1:t,[posterior.SigmaPhi]);
    
    options_online.display = out(t).options.display;

    % Define priors for parameters with past posterior
    options_online.priors = posterior(t);
    
    % TODO: remove once VBA_Iphi deal with isYout
    try
    options_online.priors.iQy = options.priors.iQy;
    catch
        options_online.priors = rmfield(options_online.priors,'iQy');
    end
    
end
