% demonstrates how to implement sparse priors with VBA
% The trick is to use a continuous transform of GLM parameters, such that
% Gaussian priors are mapped onto fat-tailed distributions. More precisely,
% the "sparse transformation" is a signed quadratic mapping (see
% sparseTransform.m). Under this transform, the priors become close to a
% Lapacian density, and the parameter VB-estimate should mimick a L1-norm
% (i.e. "lasso") estimator.
% Note: using such "sparse transformation" makes the native GLM model
% non-linear, which implies potential local minima.

% function demo_sparsePriors

close all
clear all


% 0- look at the qualitative properties of the "sparse" transformation 
mu = 0;
v = 1;
hf = figure('color',[1 1 1]);
gx = -5:1e-2:5;
pg = exp(-0.5*(gx-mu).^2/v);
pg = pg./sum(pg);
ps = exp(-abs(gx-mu)/(sqrt(v/2)));
ps = ps./sum(ps);
ha = subplot(3,2,1,'parent',hf);
plot(ha,gx,pg,'k')
title(ha,'Gaussian pdf')
ha = subplot(3,2,2,'parent',hf);
plot(ha,gx,ps,'k')
title(ha,'Laplacian pdf (sparse)')

N = 1e6;
X = mu + sqrt(v)*randn(N,1);
[n,x] = hist(X,1e2);
n = n./sum(n);
ha = subplot(3,2,3,'parent',hf);
plot(ha,x,n,'k')
set(ha,'xlim',[min(gx),max(gx)])
title(ha,'Gaussian histogram')
xlabel(ha,'x')
ylabel(ha,'p(x)')
ha = subplot(3,2,4,'parent',hf,'nextplot','add');
title(ha,'transformed Gaussian histogram')
ha2 = subplot(3,2,5,'parent',hf,'nextplot','add');
title(ha2,'sparse transformation')
ha3 = subplot(3,2,6,'parent',hf,'nextplot','add');
title(ha3,'effective regularization')
gp = [1e-1,0.25,1,10];
str = cell(length(gp),1);
col = 'brgm';
for i=1:length(gp)
    SX = sparseTransform(X,gp(i));
%     SX = normalize(SX);
    [n,x] = hist(SX,1e2);
    n = n./sum(n);
    plot(ha,x,n,col(i))
    str{i} = ['P=',num2str(gp(i))];
    sgx = sparseTransform(gx,gp(i));
    plot(ha2,gx,sgx,col(i))
    plot(ha3,x,-log(n),col(i))
end
set(ha,'xlim',[min(gx),max(gx)])
legend(ha,str)
legend(ha2,str)
xlabel(ha,'s(x)')
ylabel(ha,'p(s(x))')
xlabel(ha2,'x')
ylabel(ha2,'s(x)')
subgx = gx(1:1e2:end);
plot(ha2,subgx(subgx>=0),subgx(subgx>=0).^2,'k.')
plot(ha2,subgx(subgx<=0),-subgx(subgx<=0).^2,'k.')
plot(ha3,subgx(subgx>=0).^2,subgx(subgx>=0).^2,'k.')
plot(ha3,-subgx(subgx<=0).^2,subgx(subgx<=0).^2,'k.')

getSubplots

% 1- simulate "sparse" GLM and invert
p = 32;
n = 64;
A = randn(p,n);
phi1 = sparseTransform(X(randperm(n)),1);
sigma = 1;
y1 = A*phi1 + sqrt(sigma.^-1)*randn(p,1);

dims.n = 0;
dims.n_theta = 0;
dims.n_phi = n;
options.inG.X = A;
options.inG.sparseP = 1;
% options.priors.muPhi = 1e-2*ones(n,1);
options.checkGrads = 1;
[posterior,out] = VBA_NLStateSpaceModel(y1,[],[],@g_GLMsparse,dims,options);
set(gcf,'tag','dummy','name','"sparse" sim, "sparse" priors')
hf = figure('color',[1 1 1],'name','estimation accuracy');
ha = subplot(2,3,1,'parent',hf,'nextplot','add');
plot(ha,phi1,sparseTransform(posterior.muPhi,1),'k.')
tmp = corrcoef(phi1,sparseTransform(posterior.muPhi,1));
r = tmp(2,1);
plot(ha,[min(phi1),max(phi1)],[min(phi1),max(phi1)],'r')
title(ha,['sparse sim, sparse priors, F=',num2str(out.F),' ,r=',num2str(r)])
xlabel(ha,'simulated')
ylabel(ha,'estimated')

[posterior2,out2] = VBA_NLStateSpaceModel(y1,[],[],@g_GLM,dims,options);
set(gcf,'tag','dummy','name','"sparse" sim, Gaussian priors')
ha = subplot(2,3,2,'parent',hf,'nextplot','add');
plot(ha,phi1,posterior2.muPhi,'k.')
tmp = corrcoef(phi1,posterior2.muPhi);
r = tmp(2,1);
plot(ha,[min(phi1),max(phi1)],[min(phi1),max(phi1)],'r')
title(ha,['sparse sim, Gaussian priors, F=',num2str(out2.F),' ,r=',num2str(r)])
xlabel(ha,'simulated')
ylabel(ha,'estimated')

% dims.n_phi = dims.n_phi + 1;
% options.checkGrads = 0;
% [posterior2,out2] = VBA_NLStateSpaceModel(y1,[],[],@g_GLMsparse2,dims,options);
% set(gcf,'tag','dummy','name','"sparse" sim, Gaussian priors')
% ha = subplot(2,3,3,'parent',hf,'nextplot','add');
% plot(ha,phi1,posterior2.muPhi(1:end-1),'k.')
% tmp = corrcoef(phi1,posterior2.muPhi(1:end-1));
% r = tmp(2,1);
% plot(ha,[min(phi1),max(phi1)],[min(phi1),max(phi1)],'r')
% title(ha,['sparse sim, adaptive priors, F=',num2str(out2.F),' ,r=',num2str(r)])
% xlabel(ha,'simulated')
% ylabel(ha,'estimated')

phi_pi = pinv(A'*A)*A'*y1;
ha = subplot(2,3,3,'parent',hf,'nextplot','add');
plot(ha,phi1,phi_pi,'k.')
tmp = corrcoef(phi1,phi_pi);
r = tmp(2,1);
plot(ha,[min(phi1),max(phi1)],[min(phi1),max(phi1)],'r')
title(ha,['sparse sim, pseudo-inverse',' ,r=',num2str(r)])
xlabel(ha,'simulated')
ylabel(ha,'estimated')


% 2- simulate Gaussian GLM and invert
phi2 = X(randperm(n));
y2 = A*phi2 + sqrt(sigma.^-1)*randn(p,1);

[posterior,out] = VBA_NLStateSpaceModel(y2,[],[],@g_GLMsparse,dims,options);
set(gcf,'tag','dummy','name','Gaussian sim, "sparse" priors')
ha = subplot(2,3,4,'parent',hf,'nextplot','add');
plot(ha,phi2,sparseTransform(posterior.muPhi,1),'k.')
tmp = corrcoef(phi2,sparseTransform(posterior.muPhi,1));
r = tmp(2,1);
plot(ha,[min(phi2),max(phi2)],[min(phi2),max(phi2)],'r')
title(ha,['Gaussian sim, sparse priors, F=',num2str(out.F),' ,r=',num2str(r)])
xlabel(ha,'simulated')
ylabel(ha,'estimated')

[posterior2,out2] = VBA_NLStateSpaceModel(y2,[],[],@g_GLM,dims,options);
set(gcf,'tag','dummy','name','Gaussian sim, Gaussian priors')
ha = subplot(2,3,5,'parent',hf,'nextplot','add');
plot(ha,phi2,posterior2.muPhi,'k.')
tmp = corrcoef(phi2,posterior2.muPhi);
r = tmp(2,1);
plot(ha,[min(phi2),max(phi2)],[min(phi2),max(phi2)],'r')
title(ha,['Gaussian sim, Gaussian priors, F=',num2str(out2.F),' ,r=',num2str(r)])
xlabel(ha,'simulated')
ylabel(ha,'estimated')

phi_pi = pinv(A'*A)*A'*y2;
ha = subplot(2,3,6,'parent',hf,'nextplot','add');
plot(ha,phi2,phi_pi,'k.')
tmp = corrcoef(phi2,phi_pi);
r = tmp(2,1);
plot(ha,[min(phi2),max(phi2)],[min(phi2),max(phi2)],'r')
title(ha,['Gaussian sim, pseudo-inverse',' ,r=',num2str(r)])
xlabel(ha,'simulated')
ylabel(ha,'estimated')






