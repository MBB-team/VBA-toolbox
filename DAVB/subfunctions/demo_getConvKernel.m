% this demo estimate a convolution kernel from a dynamical system
% trajectory


clear variables
close all

% Choose basic settings for simulations
n_t = 5e2;
nu = 1;
f_fname = @f_rwl;
g_fname = @g_Id;
u = randn(nu,n_t);
alpha   = Inf;
sigma   = 1e2;
theta   = rand(nu,1);
phi     = [];
[y,x,x0,eta,e] = simulateNLSS(n_t,f_fname,g_fname,theta,phi,u,alpha,sigma,[],0);
displaySimulations(y,x,eta,e)

% return

% now estimate convolution kernel
f_fname = [];
g_fname = @g_conv0;
inG.dim.n_t = 16;
inG.dim.nu = size(u,1);
dim.n = 0;
dim.n_theta = 0;
dim.n_phi = inG.dim.n_t*inG.dim.nu +1;
options.inG = inG;
options.priors.a_sigma = sigma;
options.priors.b_sigma = 1;
% options.checkGrads = 1;
[posterior,out] = VBA_NLStateSpaceModel(y',vec(u'),f_fname,g_fname,dim,options);

plotUncertainTimeSeries(posterior.muPhi(1:end-1),sqrt(diag(posterior.SigmaPhi(1:end-1,1:end-1))))
hold on
plot(theta*(1-theta).^[0:inG.dim.n_t-1])
legend({'estimate','credible interval','theoretical'})