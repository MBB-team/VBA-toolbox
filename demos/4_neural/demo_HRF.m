
% Demo for Balloon model.
% This demo inverts a balloon model of the HRF.

close all
clear all


% Choose basic settings for simulations
TR = 2e0;                    % sampling time interval (in sec)
n_t = 60/TR;                % number of time samples (over 40 sec)
deltat = 5e-2;                 % micro-time resolution
decim = max([1,round(TR./deltat)]);
f_fname = @f_HRF2;               % Ballon model evolution function
g_fname = @g_HRF3;               % Balloon model observation function
disp('Extracting HRF params...')
% try
%     [theta,phi] = get_HRFparams(TR,deltat);
% catch
    theta = [
        0.2629
        0.1514
        -0.6285
        -0.2889
        0
        0.2646
        ];
    phi = [
        0.5770
        -0.5345
        ];
% end  
disp('Done.')
alpha = Inf;                  % simulated state noise precision
sigma = 1e2;                  % simulated data noise precision

% INPUT (in micro-time)
u = zeros(1,n_t*decim);
mT = floor([5,20]./deltat);
IU = 1e0*deltat;
u(1,mT(1):mT(2)) = IU;

% Build priors for model inversion
priors.muX0         = [0;0;0;0];
priors.SigmaX0      = 0e0*eye(4); % fix initial conditions
priors.muTheta      = 0*ones(6,1);
priors.SigmaTheta   = 1e0*eye(6,6);
priors.muPhi        = 0*ones(2,1);
priors.SigmaPhi     = 1e0*eye(2,2);
priors.a_alpha      = Inf;%1e6;
priors.b_alpha      = 0;%1e0;
priors.a_sigma      = 1e0;
priors.b_sigma      = 1e0;

% Build options and dim structures for model inversion
options.priors      = priors;
options.inF.logx2   = 0;
options.inF.deltat  = deltat;
options.inF.fullDCM = 0;
options.inF.linearized = 1;
options.inG.TE      = 0.04;
options.decim       = decim;
options.microU      = 1;
dim.n_theta         = 6;
dim.n_phi           = 2;
dim.n               = 4;

% Simulate time series of hidden states and observations
[y,x,x0,eta,e]   = VBA_simulate (...
    n_t,...
    f_fname,...
    g_fname,...
    theta,...
    phi,...
    u,...
    alpha,...
    sigma,...
    options,...
    priors.muX0);
displaySimulations(y,x,eta,e);

% Call inversion routine
[posterior,out] = VBA_NLStateSpaceModel(y,u,f_fname,g_fname,dim,options);

% Display inference results
[hres] = displayResults(posterior,out,y,x,x0,theta,phi,alpha,sigma);
set(hres,'name',['dt=',num2str(deltat),' sec (decim=',num2str(decim),')'])


% now reduce micro-time resolution to data sampling time
deltat = deltat*20;                 % micro-time resolution
decim = max([1,round(TR./deltat)]);
u = zeros(1,n_t*decim);
mT = floor([5,20]./deltat);
IU = 1e0*deltat;
u(1,mT(1):mT(2)) = IU;
options.inF.deltat  = deltat;
options.decim       = decim;
[p0,o0] = VBA_NLStateSpaceModel(y,u,f_fname,g_fname,dim,options);
set(gcf,'name',['dt=',num2str(deltat),' sec (decim=',num2str(decim),')'])

% Display inference results
[hres] = displayResults(p0,o0,y,x,x0,theta,phi,alpha,sigma);
set(hres,'name',['dt=',num2str(deltat),' sec (decim=',num2str(decim),')'])

% Re-display high-microtime VB inference
VBA_ReDisplay(posterior,out,1);
set(gcf,'name',['dt=',num2str(out.options.inF.deltat),' sec (decim=',num2str(out.options.decim),')'])

