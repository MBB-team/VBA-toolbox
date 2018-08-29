% Demo for distributed (locally multivariate) HRF


close all
clear variables


% Choose basic settings for simulations
TR          = 3;                    % sampling time interval (in sec)
n_t         = 60/TR;                % number of time samples (over 40 sec)
deltat      = 5e-2;                 % micro-time resolution
decim       = max([1,round(TR./deltat)]);
u           = zeros(1,n_t*decim);         % input
u(1,1)      = 1e0;
f_fname     = @f_HRF2;               % Ballon model evolution function
g_fname     = @g_HRF_distributed;               % Balloon model observation function
disp('Extracting HRF params...')
[theta,phi] = get_HRFparams(TR,deltat);
disp('Done.')
alpha       = Inf;                  % simulated state noise precision
sigma       = 1e6;                  % simulated data noise precision

% spatial observation parameters
inG.ind_hrf = 1:length(phi);
inG.n_reg = 1;
inG.n_phi = 4;
inG.B = randn(15,inG.n_phi);
phi = [ phi
        1
        2
        -1
        -2  ];
inG.ind_profile{1} = inG.ind_hrf(end)+1:length(phi);


% Build priors for model inversion
priors.muX0         = [0;0;0;0];
priors.SigmaX0      = 0e0*eye(4);
priors.muTheta      = 0.1*ones(length(theta),1);
priors.SigmaTheta   = 1e0*eye(length(theta));
priors.muPhi        = 0*ones(length(phi),1);
priors.SigmaPhi     = 1e0*eye(length(phi));
priors.SigmaPhi(1)  = 0;
priors.a_alpha      = Inf;%1e6;
priors.b_alpha      = 0;%1e2;
priors.a_sigma      = 1e0;
priors.b_sigma      = 1e0;

% Build options and dim structures for model inversion
options.priors      = priors;
options.inF.deltat  = deltat;
options.inF.fullDCM = 0;
options.inF.linearized = 0;
options.inG.TE      = 0.04;
options.decim       = decim;
options.microU      = 1;
options.inG         = inG;
dim.n_theta         = length(theta);
dim.n_phi           = length(phi);
dim.n               = 4;
% options.checkGrads = 1;

% Simulate time series of hidden states and observations
[y,x,x0,eta,e]   = VBA_simulate (n_t,f_fname,g_fname,theta,phi,u,alpha,sigma,options);

% Display simulated time series
displaySimulations(y,x,eta,e);
% disp('--paused--')
% pause




% Call inversion routine
% [posterior,out] = VBA_onlineWrapper(y,u,f_fname,g_fname,dim,options);
[posterior,out] = VBA_NLStateSpaceModel(y,u,f_fname,g_fname,dim,options);

% Display inference results
displayResults(posterior,out,y,x,x0,theta,phi,alpha,sigma)

% Make predictions
try
    options = out.options;
    [xs,ys,xhat,vx,yhat,vy] = VBA_comparePredictions(n_t,theta,phi,u,alpha,sigma,options,posterior,dim);
catch
    disp('------!!Unable to form predictions!!------')
end




