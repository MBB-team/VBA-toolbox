% demo for confounds removal
% This demo first simulates data from a GLM. The design matrix (X) is then
% partitionned into effects of interest and confounds (X0). Four model
% inversion are then compared:
% - full model inversion (serves as a reference)
% - reduced model inversion (without confounds)
% - reduced model inversion on pre-whitened data
% - reduced model inversion on pre-whitened data + pre-whitening of reduced
% design matrix
% Basically, one finds that the full model inversion (which is correct) is
% equivalent to the reduced model inversion on pre-whitened data +
% pre-whitening of reduced design matrix.

clear variables
close all

% Choose basic settings for simulations
ns = 16; % number of subjects
n = 256; % number of observations per subject 
g_fname = @g_GLM;
f_fname = [];
alpha   = [];
sigma   = 1e0;
theta   = [];
phi     = 8*randn(ns,1);
u       = [];

% Build options structure
X = randn(n,floor(ns/4)+1); % 1st-level design matrix
X = repmat(X,1,4);
X = X + 1e-1*randn(size(X));
X = X(:,1:ns);
inG.X = X;
options.inG     = inG;
priors.a_sigma = 1e0;
priors.b_sigma = 1e0;
priors.muPhi = zeros(ns,1);
priors.SigmaPhi = 1e4*eye(ns);
options.priors      = priors;
dim.n_theta         = 0;
dim.n_phi           = ns;
dim.n               = 0;


% Build time series of hidden states and observations
[y,x,x0,eta,e] = simulateNLSS(1,f_fname,g_fname,theta,phi,u,alpha,sigma,options,[]);


% Invert model with confounds
[p,o] = VBA_NLStateSpaceModel(y,u,f_fname,g_fname,dim,options);
displayResults(p,o,y,x,x0,theta,phi,alpha,sigma)

% remove counfounds from the model
n0 = floor(ns/2)+1:ns;
options.priors.SigmaPhi(n0,n0) = zeros(length(n0));
[p0,o0] = VBA_NLStateSpaceModel(y,u,f_fname,g_fname,dim,options);
displayResults(p0,o0,y,x,x0,theta,phi,alpha,sigma)


% pre-whiten data + remove confounds from the model
X0 = inG.X(:,n0);
% X0 = spm_orth(X0);
P0 = X0*inv(X0'*X0)*X0';
IP0 = eye(n) - P0;
y0 = IP0*y;
options.priors.SigmaPhi(n0,n0) = zeros(length(n0));
[p0,o0] = VBA_NLStateSpaceModel(y0,u,f_fname,g_fname,dim,options);
displayResults(p0,o0,y0,x,x0,theta,phi,alpha,sigma)

% pre-whiten data and observation function + remove confounds from the model
options.inG.X = IP0*inG.X;
[p1,o1] = VBA_NLStateSpaceModel(y0,u,f_fname,g_fname,dim,options);
displayResults(p1,o1,y0,x,x0,theta,phi,alpha,sigma)

