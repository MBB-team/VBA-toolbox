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


% simulation parameters
theta = [1;-2]; % volatilities
phi = log(2); % inverse temperature = 4


f_fname = @f_VBfree;
g_fname = @g_ExpUtil;
h_fname = @h_randOutcome;

% allocate feedback struture for simulations
u0 = repmat([ones(1,50),0*ones(1,50)],1,4); % 'correct' answers
fb.inH.u0 = u0;
fb.inH.er = 1; % expected reward when correct answer
fb.inH.vr = .1; % reward variance
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
priors.muTheta = [0;0];
priors.SigmaTheta = 1e0*eye(dim.n_theta);
% priors.SigmaTheta(2,2) = 0;
priors.muX0 = x0;%zeros(dim.n,1);
priors.SigmaX0 = 0e0*eye(dim.n);
priors.a_alpha = Inf;
priors.b_alpha = 0;

options.priors = priors;
options.sources = struct('type',1 ,'out', 1); % one binomial observation;
options.skipf = zeros(1,n_t);
options.skipf(1) = 1; % apply identity mapping from x0 to x1.

[y,x,x0,eta,e,u] = VBA_simulate (n_t,f_fname,g_fname,theta,phi,u,Inf,Inf,options,x0,fb);

figure
plot(y-e,'r')
hold on
plot(y,'kx')
legend({'p(y=1|theta,phi,m)','binomial data samples'})
figure
ti = {'mu1: E[1st moment of u^{o}]','s1: V[1st moment of u^{o}]','mu2: E[log- 2nd moment of u^{o}]','s2: V[log- 2nd moment of u^{o}]'};
for i=1:4
    subplot(2,2,i),plot(x([i,i+4],:)'),title(ti{i})
end
drawnow
VBA_getSubplots ();


[posterior,out] = VBA_NLStateSpaceModel(y,u,f_fname,g_fname,dim,options);

displayResults(posterior,out,y,x,x0,theta,phi,Inf,Inf);


