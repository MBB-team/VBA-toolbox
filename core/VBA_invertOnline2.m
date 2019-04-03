function [posterior, out] = VBA_invertOnline2(y, u, f_fname, g_fname, dim, options)
% This function inverts a nonlinear state-space model 'on line', i.e.
% estimate hidden states, observation/evolution parameters and
% hyperparameters in 'real-time', by adding the time samples one after the
% other. It returns a time-dependent posterior density over all unknown
% variables. See VBA_NLStateSpaceModel.m for I/O arguments.
% NB: this online wrapper does not deal with ODE state-space models.

options.OnLine = 1;

if ~ isfield(options, 'isYout')
    isYout = zeros(size(y));
else
    isYout = options.isYout;
end
%
options_online = options;

VBA_figure;
hPhi = subplot(2,1,1);
hTheta = subplot(2,1,2);


% first timestep
options_online.DisplayWin = false;
options_online.updateX0 = true;
options_online.isYout = isYout(:,1);
y_online = y(:,1);
u_online = u(:,1);

[posterior_online,out_online] = VBA_NLStateSpaceModel(y_online,u_online,f_fname,g_fname,dim,options_online);

options_online.updateX0 = false;

for t = 2 : size(y,2)
    
    t
    
    %
    % Define priors for parameters with past posterior
    options_online.priors = posterior_online(t-1);
    try
        options_online.priors.muX0 = options_online.priors.muX;
        options_online.priors.SigmaX0 = options_online.priors.SigmaX.current{1};
        options_online.priors = rmfield(options_online.priors, 'muX');
        options_online.priors = rmfield(options_online.priors, 'SigmaX');
    end
    
    % TODO: remove once VBA_Iphi deal with isYout directly
    try
    options_online.priors.iQy = options.priors.iQy(t,:);
    catch
        options_online.priors = rmfield(options_online.priors,'iQy');
    end
    try
        options_online.priors.iQx = options.priors.iQx(t);
    catch
        options_online.priors = rmfield(options_online.priors,'iQx');
    end
    %
    options_online.isYout = isYout(:,t);
    y_online = y(:,t);
    u_online = u(:,t); % TODO: microu
    
    [posterior_online(t),out_online(t)] = VBA_NLStateSpaceModel(y_online,u_online,f_fname,g_fname,dim,options_online);

    if out_online(end).options.dim.n_phi > 0
        muPhi = cat(2,out_online(1).options.priors.muPhi,[posterior_online.muPhi]);
        SigmaPhi = cat(2,diag(out_online(1).options.priors.SigmaPhi),VBA_getVar({posterior_online.SigmaPhi}));
        plotUncertainTimeSeries(muPhi,SigmaPhi,0:t,hPhi);
        xlim([0,size(y,2)])
    end
    if out_online(end).options.dim.n_theta > 0
        muTheta = cat(2,out_online(1).options.priors.muTheta,[posterior_online.muTheta]);
        SigmaTheta = cat(2,diag(out_online(1).options.priors.SigmaTheta),VBA_getVar({posterior_online.SigmaTheta}));
        plotUncertainTimeSeries(muTheta,SigmaTheta,0:t,hTheta);
        xlim([0,size(y,2)])
    end
    drawnow
   
end

posterior = posterior_online(end);
posterior.online = posterior_online;

out = out_online(end);
out.online = out_online;


VBA_ReDisplay(posterior,out);

