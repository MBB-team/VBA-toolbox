
% clear all
close all
clear variables

% class frequencies:
nk = 6; % number of classes
p = 12; % dimension of the feature space
n = 4e2; % sample size
K = 128; % max number of classes


%-- Generate n samples from a mixture of p-dimensional gaussians (MoG):

% class patterns:
mu = 1e0*randn(p,nk);
gamma= 1e-1*ones(nk,1);
lambda = ones(nk,1)./nk;
[y,labels] = generateGMM(n,lambda,mu,gamma,0);

try
    [handles] = plotResults(y,labels,mu,NaN,struct('gamma',gamma),nk,struct('verbose',1));
    set(handles.hf,'name','simulated data structure')
    drawnow
end


%-- Now classify data samples using VB inversion of GMM


% old version of the VB-MoG
[xi,eta,F,theta,K_opt] = VBEM_GM(y,K);

% display results
hf = figure('color',[1 1 1],'name','old version of VB-MoG');
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

[handles] = plotResults(y,xi,eta,F,theta,K_opt,struct('verbose',1));
set(handles.hf,'name','old version of the VB-MoG')


% new version of the VB-MoG
options.normalize = 1;
[posterior,out] = VBA_MoG(y,K,options);

[handles] = VBA_projectMoG(posterior,out,y);
set(handles.hf,'name','new version of the VB-MoG')

% display results
hf = figure('color',[1 1 1],'name','new version of VB-MoG');
h(1) = subplot(2,2,1,'parent',hf);
imagesc(mu,'parent',h(1)),colorbar
set(h(1),'xtick',1:nk,'ytick',1:p)
title(h(1),'Simulated pattern')
h(2) = subplot(2,2,2,'parent',hf);
imagesc(posterior.muEta,'parent',h(2)),colorbar
set(h(2),'xtick',1:out.dim.K,'ytick',1:p)
title(h(2),'Estimated pattern')
h(3) = subplot(2,2,3,'parent',hf);
imagesc(labels,'parent',h(3)),colorbar
set(h(3),'xtick',1:nk,'ytick',1:10:n)
title(h(3),'Simulated labels')
h(4) = subplot(2,2,4,'parent',hf);
imagesc(posterior.z','parent',h(4)),colorbar
set(h(4),'xtick',1:out.dim.K,'ytick',1:10:n)
title(h(4),'Estimated labels')



