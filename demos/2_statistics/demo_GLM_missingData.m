% demo for GLM with missing data

clear all
close all


% generate data under GLM
n = 1024; % # observations
p = 8; % # regressors
d = 16; % # missing data in GLM design matrix
s = 1e0;
X = randn(n,p);
b = randn(p,1);
g = X*b;
y = g + s*randn(size(g));
% remove data in GLM design matrix
tmp = randperm(n*p);
Xmd = X;
Xmd(tmp(1:d)) = NaN;


% set up GLM with missing data
g_fname = @g_GLM_missingData;
inG.X = Xmd;
inG.b = 1:p;
inG.md = p+1:p+d;
inG.xmd = tmp(1:d);
dim.n = 0;
dim.n_theta = 0;
dim.n_phi = p+d;
priors.muPhi = zeros(dim.n_phi,1);
priors.SigmaPhi = 1e0*eye(dim.n_phi);
options.priors = priors;
options.inG = inG;
[posterior,out] = VBA_NLStateSpaceModel(y,[],[],g_fname,dim,options);

phi = [b;VBA_vec(X(inG.xmd))];

displayResults(posterior,out,g,[],[],[],phi,1/s.^2,[]);



