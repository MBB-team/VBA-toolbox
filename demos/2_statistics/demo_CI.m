% This demo exemplifies the post-dictive interval on a simple GLM
% Let us assume that observed data are given by: Y=aX+b+e, where e are iid
% residuals. We want to portrait our uncertainty regarding "postdicted"
% data, i.e. draw something like a confidence interval around the
% regression line. 

clear all
close all

% 0- simulate data under GLM
n = 16;
X = randn(n,1);
a = 4;
b = 1;
e = 1e0*randn(n,1);
y = a*X + b + e;

% 1- fit model (here, using vaggue priors)
options.inG.X = [X,ones(n,1)];
options.priors.SigmaPhi = 1e4*eye(2);
dim.n = 0;
dim.n_phi = 2;
dim.n_theta = 0;
[post,out] = VBA_NLStateSpaceModel(y,[],[],@g_GLM,dim,options);

% 2- get Laplace postdictive density
options.priors = post;
options.priors = rmfield(options.priors,'iQy');
sx = std(X);
X0 = VBA_vec(min(X)-3*sx:1e-1:max(X)+3*sx); % extend the postdiction outside domain of fitted data
dim.p = length(X0);
dim.n_t = 1;
options.inG.X = [X0,ones(dim.p,1)];
[muy,Vy] = VBA_getLaplace([],[],@g_GLM,dim,options);

% 3- plot regression line

% Get classical CI from sampling the posterior density
N = 1e4;
sV = VBA_sqrtm (post.SigmaPhi);
phi = repmat(post.muPhi,1,N) + sV*randn(2,N);
ev = post.b_sigma./post.a_sigma;
E = 0;% sqrt(ev)*randn(length(X0),N); % add in predicted residuals? 
Y = [X0,ones(dim.p,1)]*phi + E;
%uY = prctile(Y,95,2);
uY = arrayfun(@(j) VBA_quantile(Y(j,:),.95), 1:size(Y,1))' ;
%lY = prctile(Y,5,2);
lY = arrayfun(@(j) VBA_quantile(Y(j,:),.05), 1:size(Y,1))' ;

vY = var(Y,[],2);

hf = figure('color',[1 1 1]);
ha = subplot(2,1,1,'parent',hf,'nextplot','add');
plotUncertainTimeSeries(muy',diag(Vy)',X0',ha)
plot(ha,X,y,'k.')
plot(ha,X0,uY,'r--')
plot(ha,X0,lY,'r--')
legend(ha,{'postdictive mean','posdictive STD','data','5% and 95% CI on the mean'})
xlabel(ha,'X')
ylabel(ha,'Y')
title(ha,'regression line')
ha = subplot(2,1,2,'parent',hf,'nextplot','add');
plot(ha,X0,sqrt(diag(Vy)),'k')
plot(ha,X0,sqrt(vY),'r--')
xlabel(ha,'X')
ylabel(ha,'STD[Y|y,m]')
title(ha,'postdictive standard deviation')
set(ha,'xlim',[min(X0),max(X0)])

VBA_getSubplots ();
