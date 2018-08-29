% demo for log-normal priors

clear all
close all


g_fname = @g_exp;  % gx = phi(2)*exp(phi(1)*in.x)+phi(3);
f_fname = [];
priors.muPhi = [0;1;0];
priors.SigmaPhi = zeros(3,3);
priors.SigmaPhi(1,1) = 0.1;
priors.a_alpha = 0;
priors.b_alpha = 0;
priors.a_sigma = Inf;
priors.b_sigma = 0;
inG.x = 1;

u = [];
options.inG = inG;
options.priors = priors;
n_t = 1;
dim.n = 0;
dim.n_theta = 0;
dim.n_phi = 3;
dim.p = 1;
N = 2^11;

np = 64;
lx = [];
ly = [0,1e1];
[pX,gX,pY,gY,X,Y] = VBA_MCMC_predictiveDensity(f_fname,g_fname,u,n_t,options,dim,N,np,lx,ly);


figure('color',[1 1 1])
plot(gY,pY./sum(pY),'marker','.')
title('log-normal probability density function')
xlabel('x')
ylabel('p(x)')