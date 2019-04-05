function [posterior, out] = demo_multisession ()
% // VBA toolbox //////////////////////////////////////////////////////////
%
% [posterior, out] = demo_multisession ()
% Shows hot to use to fit multiple sessions at once
%
% /////////////////////////////////////////////////////////////////////////

%% Model definition 
% =========================================================================

% evolution function
function [fx] = f_demo_multisession (Xt, Theta, ut, ~)
    % x(t+1) = x(t) + u(t) theta(1) + theta(2)
    fx = Xt + [ut 1] * Theta;
end

function [gx] = g_demo_multisession (Xt, Phi, ~, ~)
    % y(t) = x(t) + phi(1) + noise 
    gx = Phi + Xt;
end


%% Simulations
% =========================================================================
% To emulate multisession, we simulate the same model twice with different 
% set of parameters

% total number of observation
n_t = 100;
% inputs for overal experiment
u = rand (1, n_t);
% no display
options.verbose = false;

% first session
theta = [5; 1];
phi = 5;
X0 = 0;
y1 = VBA_simulate (round(n_t/2),@f_demo_multisession,@g_demo_multisession,theta,phi,u(1:(n_t/2)),Inf,.01,options,X0);

% second session
theta = - theta; % reverse evolution
phi = phi - 3; % shift in observations
X0 = 0;
y2 = VBA_simulate (round(n_t/2),@f_demo_multisession,@g_demo_multisession,theta,phi,u((1+n_t/2):end),Inf,.01,options,X0);

% concatenate both sessions
y = [y1 y2];

%% Model inversion
% =========================================================================

% dimension of the problem
dim.n_t = n_t; % number of observation timepoints
dim.n_phi = 1; % nb of observation parameters
dim.n_theta = 2; % nb of evolution parameters
dim.n = 1; % nb of states

% display options
options.verbose = true;

% fit data as 2 sessions
% -------------------------------------------------------------------------
% split into sessions
options.multisession.split = [n_t/2 n_t/2]; % two sessions of 120 datapoints each

% By default, all parameters are duplicated for each session. However, you
% can fix some parameters so they remain constants across sessions.

% + Example: same evolution parameter in both sessions
% options.multisession.fixed.theta = 1; % <~ uncomment for fixing theta(1)

% + Example: same observation parameter in both sessions
% options.multisession.fixed.phi = 1; % <~ uncomment for fixing phi(1)

% + Example: same initial state in both sessions
% options.multisession.fixed.X0 = 1; % <~ uncomment for fixing X0(1)

% Model identification as usual
[posterior.split,out.split] = VBA_NLStateSpaceModel(y,u,@f_demo_multisession,@g_demo_multisession,dim,options);
set(out.split.options.hf, 'name','demo_multisession: duplicated parameters')

% 1 session case
% -------------------------------------------------------------------------
% fix all parameters to equal across sessions. Note that it would be
% simpler to just remove the multisession options!
options.multisession.fixed.theta = 'all'; % shrotcut for 1:dim.n_theta
options.multisession.fixed.phi   = 'all'; % shrotcut for 1:dim.n_phi

% Model identification as usual
[posterior.fixed,out.fixed] = VBA_NLStateSpaceModel(y,u,@f_demo_multisession,@g_demo_multisession,dim,options);
set(out.fixed.options.hf, 'name','demo_multisession: fixed parameters')

end








