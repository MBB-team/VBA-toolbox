function [posterior_split,out_split,posterior_fixed,out_fixed] = demo_multisession()


%% ### Model Definition 

function [fx]=f_demo_multisession(Xt,Theta,ut,~)
fx = Xt + [ut 1]*Theta;
% dfdx = 1;
% dfdp = ut;
end

function [gx]=g_demo_multisession(Xt,Phi,~,~)
gx = Phi+Xt;
% dgdx = 1;
% dgdp = 1;
end


%% ### Simulations

% We simulate the same model two times with different parameters
n_t = 100;
u = rand(1,n_t);
options.verbose=0;

% first session
theta = [5 1]';
phi = 5;
X0 = 0;
y1 = simulateNLSS(round(n_t/2),@f_demo_multisession,@g_demo_multisession,theta,phi,u(1:(n_t/2)),Inf,.1,options,X0);

% second session
theta = -[5 1]';
phi = 2;
X0 = 0;
y2 = simulateNLSS(round(n_t/2),@f_demo_multisession,@g_demo_multisession,theta,phi,u((1+n_t/2):end),Inf,.1,options,X0);

% concatenate both sessions
y = [y1 y2];

%% ### Model inversion

% = Set the priors for as if data were from a unique session
priors.muX0       = 0;
priors.SigmaX0    = .1;
priors.muTheta    = [0 0]';
priors.SigmaTheta = eye(2);     
priors.muPhi      = 0;       % prior mean on observation params
priors.SigmaPhi   = 1;       % prior covariance on observation params
priors.a_sigma    = 1;       % Jeffrey's prior
priors.b_sigma    = 1;       % Jeffrey's prior   
options.priors    = priors;  % include priors in options structure

dim.n_t = n_t;
dim.n_phi = 1;            % nb of observation parameters
dim.n_theta = 2;          % nb of evolution parameters
dim.n=1;  


options.verbose=1;
options.DisplayWin=1;
options.extended=1
% =========================================================================
% 2 sessions case
% =========================================================================

% split into sessions
options.multisession.split = [n_t/2 n_t/2]; % two sessions of 120 datapoints each

% % By default, all parameters are duplicated for each session. However, you
% % can fix some parameters so they remain constants across sessions.

% + Example: same evolution parameter in both sessions
% options.multisession.fixed.theta = 1; % <~ uncomment for fixing theta(1)

% + Example: same observation parameter in both sessions
% options.multisession.fixed.phi = 1; % <~ uncomment for fixing phi(1)

% + Example: same initial state in both sessions
% options.multisession.fixed.X0 = 1; % <~ uncomment for fixing X0(1)

% = Model identification as usual
[posterior_split,out_split] = VBA_NLStateSpaceModel(y,u,@f_demo_multisession,@g_demo_multisession,dim,options);

pause()

% =========================================================================
% 1 session case
% =========================================================================

% keep the two sessions but fix all parameters
options.multisession.fixed.theta = 'all'; % shrotcut for 1:dim.n_theta
options.multisession.fixed.phi   = 'all'; % shrotcut for 1:dim.n_phi

% = Model identification as usual
[posterior_fixed,out_fixed] = VBA_NLStateSpaceModel(y,u,@f_demo_multisession,@g_demo_multisession,dim,options);

end








