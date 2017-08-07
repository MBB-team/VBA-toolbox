% This script simulates and inverts a Bayesian sequence learner (BSL)

clear all
close all
clc

% simulation parameters
K = 1;
options.inF = struct('K',K);
options.inG = struct('K',K);
dim.n = 2^(K+1);
dim.n_theta = 1;
dim.n_phi = 2;



%% simulate sequence of BSL choices
x0 = 1*ones(dim.n,1); % log-odds of P(o=1)
theta = [-2]; % BSL's prior volatility
phi = [-1;0]; % (log-) inverse-temperature, bias
N = 150; % number of trials
p = 0.75;
P = repmat([p,1-p,p,1-p],1,N); % probabilistic repetition of [1 0 0 1]
% P = repmat(p,1,N); % probabilistic repetition of [1 0 0 1]
y = zeros(1,N);
for i=1:N
    tmp = VBA_sample('multinomial',struct('p',[P(i);1-P(i)],'n',1),1,0);
    y(i) = tmp(1);
end
a = zeros(1,N);
gx = NaN(1,N);
x = zeros(dim.n,N);
x(:,1:K) = repmat(x0,1,K); %initialize hidden states
u(:,1:K) = NaN(K+1,K);
for i=K+2:N
    u(:,i) = flipud(vec(y(i-K-1:i-1)));
    if K==0 && i==1 % the issue only arises for 0-BSL (degenerated!)
        x(:,i) = f_BSL(x0,theta,u(:,i),options.inF);
    else
        x(:,i) = f_BSL(x(:,i-1),theta,u(:,i),options.inF);
    end
    gx(i) = g_BSL(x(:,i),phi,u(:,i),options.inG);
    tmp = VBA_sample('multinomial',struct('p',[gx(i);1-gx(i)],'n',1),1,0);
    a(i) = tmp(1);%gx(i)>.5;
end
hf = figure('color',[1 1 1]);
ha = subplot(2,1,1,'parent',hf);
plot(ha,x')
ha = subplot(2,1,2,'parent',hf);
plot(ha,gx)
hold on, plot([1,N],[0.5,0.5],'color',0.2*[1 1 1])
ic = find(a(1:N)==y(1:N));
plot(ha,ic,a(ic),'g.')
ii = find(a(1:N)~=y(1:N));
plot(ha,ii,a(ii),'r.')
title(ha,['perf=',num2str(100*length(ic)./N,3),'%'])


%% invert BSL model given sequence of agent's choices
options.skipf = zeros(1,N);
options.skipf(1:K+1) = 1;
options.sources = struct('type',1 ,'out', 1); % one binomial observation;
options.SigmaTheta = 1e2*eye(dim.n_theta);
f_fname = @f_BSL;
g_fname = @g_BSL;
[posterior,out] = VBA_NLStateSpaceModel(a,u,f_fname,g_fname,dim,options);

%% Display results
displayResults(posterior,out,a,x,x0,theta,phi,[],[])
hf = unwrapKBSL(posterior.muX,posterior.muPhi,u,options.inG);

