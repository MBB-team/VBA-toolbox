% Demo for DCM for fMRI (without Balloon observer)
% This demo inverts the DCM for fMRI model, without the ballon model, which
% is replaced by a nonlinear sigmoid observation function.


% close all
clear variables



% Choose basic settings for simulations
n_t = 1.5e2;                      % number of time samples
TR = 1e0;                       % sampling period (in sec)
microDT = 1e-1;                 % micro-time resolution (in sec)
nreg = 3;
nu = 2;
u       = zeros(2,n_t);         % inputs
u(1,2:3) = 1;
u(1,30:31) = 1;
u(1,50:51) = 1;
u(2,30:31) = 1;
u(2,60:61) = 1;
% u       = zeros(2,n_t);
% u(1,100:150) = 1;
% u(1,300:350) = 1;
% % u(1,500:700) = 1;
% u(2,300:350) = 1;
% % u(2,600:650) = 1;
alpha   = 1e2;
sigma   = 1e3;
theta   = [ -0.5
            -0.5
            -0.5
            -2.5
            -1.5
            -.5
             .1
            -2
            -0.1];
theta = exp(theta);
theta(end) = log(theta(end));
phi     = [];

% DCM specification
A = [0 1 1
     1 0 1
     0 1 0];
B{1} = zeros(nreg,nreg);
B{2} = [0 0 0
        1 0 0
        0 0 0];
C = [1 0
     0 0
     0 0];
D{1} = [0 0 0
        0 0 0
        0 1 0];
D{2} = zeros(nreg,nreg);
D{3} = zeros(nreg,nreg);


% Prepare precalculated matrices for DCM inversion
[inF] = prepare_dcm(A,B,C,D);
inG.G0 = 50;
inG.k = 1;

% Build priors for model inversion
priors.muX0 = 0*ones(nreg,1);
priors.SigmaX0 = 0e-1*eye(nreg);
priors.muTheta = 0*ones(size(theta));
priors.muTheta(inF.indself) = -0;
priors.SigmaTheta = 1e-1*eye(length(theta));
priors.a_alpha = 1e0;
priors.b_alpha = 1e0;
priors.a_sigma = 1e0;
priors.b_sigma = 1e0;

% Build options and dim structures for model inversion
priors.iQx = cell(n_t,1);
for t=1:n_t
    priors.iQx{t} = eye(nreg);
    priors.iQx{t}(3,3) = 1e4;
end

% Build options structure for temporal integration of SDE
options.priors      = priors;
options.inF     = inF;
options.inG     = inG;
options.decim = max([1,floor(TR./microDT)]);
options.inF.deltat = TR./options.decim;
options.GnFigs = 1;
options.Laplace = 0;
options.gradF = 0;

dim.n_theta         = length(theta);
dim.n_phi           = 0;
dim.n               = nreg;

f_fname = @f_dcm4fmri;
g_fname = @g_Id;

% Build time series of hidden states and observations
[y,x,x0,eta,e] = simulateNLSS(n_t,f_fname,g_fname,theta,phi,u,alpha,sigma,options);

% display time series of hidden states and observations
displaySimulations(y,x,eta,e)
% disp('--paused--')
% pause


% Call inversion routine
% u = zeros(size(u));
% options.embed = 1;
[posterior,out] = VBA_NLStateSpaceModel(y,u,f_fname,g_fname,dim,options);

% Display results
displayResults(posterior,out,y,x,x0,theta,phi,alpha,sigma)

% Make predictions
try
    options = out.options;
    [xs,ys,xhat,vx,yhat,vy] = comparePredictions(...
        n_t,theta,phi,u,alpha,sigma,options,posterior,dim);
catch
    disp('------!!Unable to form predictions!!------')
end

