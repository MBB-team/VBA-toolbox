% demo operant learning with VB volatile observer (adapted from Mathys, Daunizeau et al. 2010)


close all
clear variables
clc



% evolution, observation and feedback functions
f_fname = @f_OpLearn;
g_fname = @g_VBvolatile0;
h_fname = @h_truefalse;


% allocate feedback struture for simulations
u0 = [ones(1,50)]; % possible feedbacks
fb.inH.u0 = repmat([u0,~u0],1,6); % with reversals
fb.h_fname = h_fname;
fb.indy = 1;
fb.indfb = 2;


% simulation parameters
theta = [1;-4;-1];
inF.lev2 = 0; % remove 3rd level (volatility learning)
inF.kaub = 1.4;
inF.thub = 1;
inF.rf = -1;
phi = [1;0.5]; % inverse temperature = 2, bias towards 'a=1' = 1
inG.respmod = 'taylor';

% choose initial conditions
x0 = repmat([0.5;0;0;1;log(4)],2,1);
u = zeros(2,size(fb.inH.u0,2)+1);

dim = struct('n',2*5,'n_theta',3,'n_phi',2);

priors.muPhi = zeros(dim.n_phi,1);
priors.muTheta = [0;-4;0];
priors.muX0 = x0;
priors.SigmaPhi = 1e2*eye(dim.n_phi);
priors.SigmaTheta = 1e2*eye(dim.n_theta);
priors.SigmaX0 = 0e1*eye(dim.n);
priors.a_alpha = Inf;
priors.b_alpha = 0;

options.priors = priors;
options.binomial = 1;
options.inF = inF;
options.inG = inG;
options.skipf = zeros(1,length(u));
options.skipf(1) = 1; % apply identity mapping from x0 to x1.

[y,x,x0,eta,e,u] = simulateNLSS_fb(length(u),f_fname,g_fname,theta,phi,u,Inf,Inf,options,x0,fb);

figure
plot(y-e,'r')
hold on
plot(y,'kx')
legend({'p(y=1|theta,phi,m)','binomial data samples'})
getSubplots
pause


% options.isYout = zeros(1,size(y,2));
% options.isYout(75:125) = 1;

[posterior,out] = VBA_NLStateSpaceModel(y,u,f_fname,g_fname,dim,options);

displayResults(posterior,out,y,x,x0,theta,phi,Inf,Inf)


