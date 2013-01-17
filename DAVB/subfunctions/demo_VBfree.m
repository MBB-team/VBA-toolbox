% demo VB-learning (free learner)
% This demo simulates a VB agent that reacts to outcomes by choosing the
% action that maximizes the expected gain. Here, the agent learns the
% distribution of outcomes associated with each available action. The key
% trick is that this distribution is summmarized in terms of its two first
% moments. The agent a priori belives that these can drift over time, with
% transition variances (i.e. volatilities) exp(theta), where theta is the
% 2x1 evolution parameter vector.
% NB: there are 2 posteriro moments per moment of the outcome distribution,
% which is action-dependent. With 2 available actions, this means there are
% 2x2x2=8 hidden states in this model (to be compared with 2 Q-values for
% RL).
% NB2: the action emission law does not use any utility mapping yet!


close all
clear variables
clc


% simulation parameters
theta = [1;-32]; % volatilities
phi = log(4); % inverse temperature = 4


f_fname = @f_VBfree;
g_fname = @g_ExpUtil;
h_fname = @h_truefalse;

% allocate feedback struture for simulations
u0 = [ones(1,50)]; % possible feedbacks
fb.inH.u0 = [u0,~u0,u0,~u0,u0,~u0,u0]; % with reversals
fb.h_fname = h_fname;
fb.indy = 1;
fb.indfb = 2;

% choose dummy initial conditions
x0 = repmat([0;1;0;1],2,1);
u = zeros(2,size(fb.inH.u0,2)+1);

n_t = length(u); % number of trials

dim = struct('n',8,'n_theta',2,'n_phi',1);

priors.muPhi = zeros(dim.n_phi,1);
priors.SigmaPhi = 1e0*eye(dim.n_phi);
priors.muTheta = [0;-32];%0*ones(dim.n_theta,1);
% priors.muTheta(2) = theta(2);
priors.SigmaTheta = 1e0*eye(dim.n_theta);
priors.SigmaTheta(2,2) = 0;
priors.muX0 = zeros(dim.n,1);
priors.SigmaX0 = 1e0*eye(dim.n);
priors.a_alpha = Inf;
priors.b_alpha = 0;

options.priors = priors;
options.binomial = 1;
options.skipf = zeros(1,n_t);
options.skipf(1) = 1; % apply identity mapping from x0 to x1.

[y,x,x0,eta,e,u] = simulateNLSS_fb(n_t,f_fname,g_fname,theta,phi,u,Inf,Inf,options,x0,fb);

figure
plot(y-e,'r')
hold on
plot(y,'kx')
legend({'p(y=1|theta,phi,m)','binomial data samples'})
figure
ti = {'mu1','s1','mu2','s2'};
for i=1:4
    subplot(2,2,i),plot(x([i,i+4],:)'),title(ti{i})
end
drawnow
getSubplots

% pause


% options.isYout = zeros(1,size(y,2));
% options.isYout(75:125) = 1;

[posterior,out] = VBA_NLStateSpaceModel(y,u,f_fname,g_fname,dim,options);

displayResults(posterior,out,y,x,x0,theta,phi,Inf,Inf)


