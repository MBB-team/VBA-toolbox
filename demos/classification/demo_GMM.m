
% clear all
close all
clear variables

% class frequencies:
nk = 3; % number of classes
p = 2; % dimension of the feature space
n = 4e2; % sample size
K = 32;ceil(n/2); % max number of classes


%-- Generate n samples from a mixture of p-dimensional gaussians (MoG)
mu = 1e3+ 1e3*randn(p,nk); % class means
gamma= 1e5*ones(nk,1); % class variances
lambda = ones(nk,1)./nk; % class frequencies
[y,labels] = generateGMM(n,lambda,mu,gamma,0);

try
    [handles] = plotResults(y,labels,mu,NaN,struct('gamma',gamma),nk,struct('verbose',1));
    set(handles.hf,'name','simulated data structure')
    drawnow
end


%-- Now classify data samples using VB inversion of GMM

% NB: data range matters when performing VB clustering. The following
% specifications does not make VB robust for any data structure!
options.normalize = 1;
options.MaxIter = 128;
for k=1:K
    options.priors.a_gamma = 1e0*ones(K,1);
    options.priors.b_gamma = 1e-2*ones(K,1);
    options.priors.SigmaEta{k} = p*1e-2*eye(p);
end
[posterior,out] = VBA_MoG(y,K,options);

%-- Display results

[handles] = VBA_projectMoG(posterior,out,y);
try,set(handles.hf,'name','new version of the VB-MoG'),end

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



