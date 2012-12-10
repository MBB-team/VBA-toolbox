
% clear all
close all
clear variables

% 
% path_sampling = 'C:\Users\JeanD\Documents\MATLAB\classification\sampling';%[fileparts(mfilename('fullpath')),filesep,'..',filesep,'sampling']; %'C:\Users\jdaunizeau\Documents\MatlabWork\Routinetheque\sampling';
% addpath(path_sampling)

% class frequencies:
nk = 12; % number of classes
p = 8; % dimension of the feature space
n = 5e2; % sample size
lambda = ones(nk,1);
lambda = lambda./sum(lambda);

% class patterns:
mu = 1e0*randn(p,nk);
% mu(:,end) = 1e1*randn(p,1);
gamma= 1e-1*ones(nk,1);

% Generate n data samples from a binomial mixture model (BMM):
[y,labels] = generateGMM(n,lambda,mu,gamma,0);




% classify data samples using VB inversion of BMM
K = 32; % max number of classes


[xi,eta,F,theta,K_opt] = VBEM_GM(y,K)


[out.handles] = plotResults(y,xi,eta,F,theta,K_opt,struct('verbose',1));


options.verbose = 1;
% options.priors.muEta = mu;
% options.init = 'prior';
options.minSumZ = 1e-2;
% options.priors.d = ones(K,1)./K;
[posterior,out] = VBA_MoG(y,K,options);


xi = posterior.z';
eta = posterior.muEta;
F = out.F(end);
K_opt = size(posterior.z,1);


[out.handles] = plotResults(y,xi,eta,F,struct('gamma',posterior.a_gamma./posterior.b_gamma),K_opt,struct('verbose',1));

% display results
hf = figure('color',[1 1 1]);
h(1) = subplot(2,2,1,'parent',hf);
imagesc(mu,'parent',h(1)),colorbar
set(h(1),'xtick',1:nk,'ytick',1:p)
title(h(1),'Simulated pattern')
h(2) = subplot(2,2,2,'parent',hf);
imagesc(eta,'parent',h(2)),colorbar
set(h(2),'xtick',1:K_opt,'ytick',1:p)
title(h(2),'Estimated pattern')
h(3) = subplot(2,2,3,'parent',hf);
imagesc(labels,'parent',h(3)),colorbar
set(h(3),'xtick',1:nk,'ytick',1:10:n)
title(h(3),'Simulated labels')
h(4) = subplot(2,2,4,'parent',hf);
imagesc(xi,'parent',h(4)),colorbar
set(h(4),'xtick',1:K_opt,'ytick',1:10:n)
title(h(4),'Estimated labels')



