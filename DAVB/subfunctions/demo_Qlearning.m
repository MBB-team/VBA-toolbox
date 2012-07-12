% demo Q-learning


close all
clear variables
clc


% simulation parameters
theta = sigm(0.3,struct('INV',1)); % learning rate = 0.3
phi = log(2); % inverse temperature = 2


f_fname = @f_Qlearn2;
g_fname = @g_softmax;
h_fname = @h_truefalse;

% allocate feedback struture for simulations
u0 = [ones(1,100)]; % possible feedbacks
fb.inH.u0 = [u0,~u0,u0,~u0,u0,~u0]; % with reversals
fb.h_fname = h_fname;
fb.indy = 1;
fb.indfb = 2;

% choose dummy initial conditions
x0 = zeros(2,1);
u = zeros(2,size(fb.inH.u0,2)+1);

n_t = length(u); % number of trials

dim = struct('n',2,'n_theta',1,'n_phi',1);

priors.muPhi = zeros(dim.n_phi,1);
priors.muTheta = ones(dim.n_theta,1);
priors.muX0 = zeros(2,1);
priors.SigmaPhi = 1e0*eye(dim.n_phi);
priors.SigmaTheta = 1e0*eye(dim.n_theta);
priors.SigmaX0 = 0e1*eye(dim.n);
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
getSubplots
% pause


% options.isYout = zeros(1,size(y,2));
% options.isYout(75:125) = 1;

[posterior,out] = VBA_NLStateSpaceModel(y,u,f_fname,g_fname,dim,options);

displayResults(posterior,out,y,x,x0,theta,phi,Inf,Inf)


