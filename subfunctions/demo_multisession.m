function [posterior1,out1,posterior2,out2,posterior,out,posterior_fixed,out_fixed] = demo_multisession()


%% ### Model Definition 

function [fx]=f_demo_multisession(Xt,Theta,ut,inF)
fx = Xt + [ut 1]*Theta;
% dfdx = 1;
% dfdp = ut;
end

function [gx]=g_demo_multisession(Xt,Phi,ut,inG)
gx = Phi+Xt;
% dgdx = 1;
% dgdp = 1;
end


%% ### Simulations
% We simulate the same model two times with different parameters
n_t = 150;
u = 3*rand(1,n_t);
options.verbose=0;

% first session
theta = [5 1]';
phi = 5;
X0 = 0;
y1 = simulateNLSS(round(n_t/2),@f_demo_multisession,@g_demo_multisession,theta,phi,u(1:(n_t/2)),Inf,.1,options,X0);

% second session
theta = [5 1]';
phi = 5;
X0 = 0;
y2 = simulateNLSS(round(n_t/2),@f_demo_multisession,@g_demo_multisession,theta,phi,u((1+n_t/2):end),Inf,.1,options,X0);

% concatenate both sessions
y = [y1 y2];

%% ### Model inversion

% = Set the priors for as if data were from a unique session
priors.muX0 = 0;
priors.SigmaX0 = 0;
priors.muTheta = [0 0]';
priors.SigmaTheta = eye(2);     
priors.muPhi = 0;         % prior mean on observation params
priors.SigmaPhi = 1;      % prior covariance on observation params
priors.a_sigma = 1;       % Jeffrey's prior
priors.b_sigma = 1;       % Jeffrey's prior
priors.a_alpha = Inf;     
priors.b_alpha = 0;       
options.priors = priors;  % include priors in options structure

dim.n_t = n_t;
dim.n_phi = 1;            % nb of observation parameters
dim.n_theta = 2;          % nb of evolution parameters
dim.n=1;  


options.verbose=1;
options.DisplayWin=0;

% options.TolFun = 1e-12;
% options.GnTolFun = 1e-12;

dim.n_t = n_t/2;

[posterior1,out1] = VBA_NLStateSpaceModel(y1,u(1:(n_t/2)),@f_demo_multisession,@g_demo_multisession,dim,options);
options2 = options;
options2.priors.a_sigma = posterior1.a_sigma;
options2.priors.b_sigma = posterior1.b_sigma;
[posterior2,out2] = VBA_NLStateSpaceModel(y2,u((1+n_t/2):end),@f_demo_multisession,@g_demo_multisession,dim,options2);

dim.n_t = n_t;

% = Specifiy the sessions
options.multisession.split = [n_t/2 n_t/2]; % two sessions of 120 datapoints each
% % By default, all parameters are duplicated for each session. However, you
% % can fix some parameters so they remain constants across sessions.

% + Example: same evolution parameter in both sessions
% options.multisession.fixed.theta = 1:2; % <~ uncomment for fixing theta(1)

% + Example: same observation parameter in both sessions
% options.multisession.fixed.phi = 1; % <~ uncomment for fixing phi(1)

% + Example: same initial state in both sessions
% options.multisession.fixed.X0 = 1; % <~ uncomment for fixing X0(1)

% = Model identification as usual
[posterior,out] = VBA_NLStateSpaceModel(y,u,@f_demo_multisession,@g_demo_multisession,dim,options);

options.multisession.fixed.theta = 1:2; % <~ uncomment for fixing theta(1)
options.multisession.fixed.phi = 1; % <~ uncomment for fixing phi(1)
options.multisession.fixed.X0 = 1; % <~ uncomment for fixing X0(1)

% = Model identification as usual
fprintf('-----\n');
[posterior_fixed,out_fixed] = VBA_NLStateSpaceModel(y,u,@f_demo_multisession,@g_demo_multisession,dim,options);




% = Display session per session parameter estimates
% for i=1:2
%     fprintf('=== Session %d ================== \n',i);
%     disp(posterior.perSession(i));
%     fprintf('\n');
% end

end








