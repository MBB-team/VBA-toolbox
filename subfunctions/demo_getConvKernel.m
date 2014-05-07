% this script demonstrates how to estimate a convolution kernel from a
% dynamical system trajectory. In this particular example, the dynamical
% system is in fact a Rescorla-Wagner learning agent.


clear variables
close all

% Choose basic settings for simulations
n_t = 1e3;
nu = 1;
f_fname = @f_rwl;
g_fname = @g_Id;
u = (2*randn(nu,n_t)>0)-1;%randn(nu,n_t);%
alpha   = 1e2;
sigma   = 1e2;
theta   = 0.2;%rand(nu,1);
phi     = [];


fb.inH.u = u; % with reversals
fb.h_fname = @h_Id;
fb.indy = [];
fb.indfb = 1;
options.skipf = zeros(1,n_t);
options.skipf(1) = 1; % apply identity mapping from x0 to x1.

[y,x,x0,eta,e,u] = simulateNLSS_fb(n_t,f_fname,g_fname,theta,phi,u,alpha,sigma,options,0,fb);
displaySimulations(y,x,eta,e)

dim.n = 1;
dim.n_theta = 1;
dim.n_phi = 0;
% options.priors.a_alpha = Inf;
% options.priors.b_alpha = 0;
options.backwardLag = 16;
[p0,o0] = VBA_NLStateSpaceModel(y,u,f_fname,g_fname,dim,options);


% now estimate convolution kernel
g_fname = @g_conv0;
inG.dim.n_t = 16;
inG.dim.nu = size(u,1);
dim.n = 0;
dim.n_theta = 0;
dim.n_phi = inG.dim.n_t*inG.dim.nu +1;
opt.inG = inG;
% opt.priors.a_sigma = sigma;
% opt.priors.b_sigma = 1;
% options.checkGrads = 1;
[posterior,out] = VBA_NLStateSpaceModel(y',vec(u'),[],g_fname,dim,opt);

plotUncertainTimeSeries(posterior.muPhi(1:end-1),sqrt(diag(posterior.SigmaPhi(1:end-1,1:end-1))))
hold on
plot(theta*(1-theta).^[0:inG.dim.n_t-1])
legend({'estimate','credible interval','theoretical'})

VBA_ReDisplay(p0,o0,1)
