function [] = demo_dcm4fmri_old()

% [broken] Demo for DCM for fMRI (with HRF embedding)
% This demo inverts the DCM for fMRI model, which contains the ballon model
% as a generalized observation function (not affected by stochastic
% innovations).

warning on
warning('This demo is now obsolete. Use demo_dcm4fmri.m!')
return

close all
clear variables

% Choose basic settings for simulations
n_t = 1e2;
f_fname = @f_fullDCM4fmri;
g_fname = @g_fullDCM4fmri;
nu = 2;
decim = 20;
deltat = 10e-2;
u       = zeros(2,n_t);
u(1,1:5) = 1;
u(1,30:35) = 1;
u(1,50:70) = 1;
u(2,30:35) = 1;
u(2,60:65) = 1;
alpha   = Inf;%1e8;
sigma   = 1e8;

% DCM specification
A = [0 1 1
     1 0 1
     0 1 0];
nreg = size(A,1);
B{1} = zeros(nreg,nreg);
B{2} = [0 0 0
        1 0 0
        0 0 0];
C = [1 0
     0 0
     0 0];
D{1} = zeros(nreg,nreg);
% D{1} = [0 0 0
%         0 0 0
%         0 1 0];
D{2} = zeros(nreg,nreg);
D{3} = zeros(nreg,nreg);

theta   = [exp([ -1.5
            -1.5
            -0.5
            -2.5
            -1.5
            -.5
             .1]);
            -1];
phi     = zeros(6*nreg,1);
phi(2)  = 0;

% Build priors for model inversion
priors.muX0 = 0*ones(nreg,1);
priors.SigmaX0 = 0*speye(nreg);
priors.muTheta = 0*ones(size(theta));
priors.SigmaTheta = 1e-1*speye(length(theta));
priors.muPhi = 0*ones(size(phi));
priors.SigmaPhi = 1e-1*speye(length(phi));
priors.a_alpha = 1e0;
priors.b_alpha = 1e0;
priors.a_sigma = 1e4;
priors.b_sigma = 1e0;


% Prepare precalculated matrices for DCM inversion
[options] = prepare_fullDCM_old(A,B,C,D,deltat,n_t,decim,priors);
dim.n_theta         = length(theta);
dim.n_phi           = length(phi);
dim.n               = nreg;



% Build time series of hidden states and observations
[y,x,x0,eta,e] = simulateNLSS(n_t,f_fname,g_fname,theta,phi,u,alpha,sigma,options);

% display time series of hidden states and observations
displaySimulations(y,x,eta,e)
disp('--paused--')
pause



% Call inversion routine
[posterior,out] = VBA_NLStateSpaceModel(y,u,f_fname,g_fname,dim,options);

% Display results
displayResults(posterior,out,y,x,x0,theta,phi,alpha,sigma)

% Make predictions
try
    options = out.options;
    [x,y,xhat,vx,yhat,vy] = comparePredictions(...
        1e3,theta,phi,zeros(size(u,1),1e3),alpha,sigma,options,posterior,dim);
catch
    disp('------!!Unable to form predictions!!------')
end


