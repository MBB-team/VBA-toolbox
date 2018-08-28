% demo operant learning with VB volatile observer (adapted from Mathys,
% Daunizeau et al. 2010)
% This demo simulates a "volatile learner" experiencing feedback to her
% chosen actions in a 2-armed bandit task. VBA estimation then proceeds.

close all
clear variables
clc



% evolution, observation and feedback functions
f_fname = @f_OpLearn;
g_fname = @g_VBvolatile0;
h_fname = @h_truefalse;


% allocate feedback struture for simulations
u0 = [randn(1,50)>-0.5]; % possible feedbacks
fb.inH.u0 = repmat([u0,1-u0],1,4); % with reversals
fb.h_fname = h_fname;
fb.indy = 1;
fb.indfb = 2;


% simulate VB volatile learner in a 2-armed bandit task
theta = [0;-2;0]; % [1;-3;-1];
phi = [0;0]; % inverse (log-) temperature & bias
inF.lev2 = 1; % 3rd level (volatility learning)
inF.kaub = 1.4;
inF.thub = 1;
inF.rf = 1;
inG.respmod = 'fixedForm'; % or 'taylor';
x0 = repmat([0;0;0;0;0],2,1);
u = zeros(2,size(fb.inH.u0,2)+1);
options.sources = struct('type',1 ,'out', 1); % one binomial observation;
options.inF = inF;
options.inG = inG;
options.skipf = zeros(1,length(u));
options.skipf(1) = 1; % apply identity mapping from x0 to x1.
[y,x,x0,eta,e,u] = VBA_simulate (length(u),f_fname,g_fname,theta,phi,u,Inf,Inf,options,x0,fb);

figure
plot(y-e,'r')
hold on
plot(y,'kx')
legend({'p(y=1|theta,phi,m)','binomial data samples'})
VBA_getSubplots ();
% pause

dummy.options = options;
[ha,hf] = unwrapVBvolatileOTO(struct('muX',x,'muTheta',theta),dummy);
set(hf,'name','simulated volatile VB learner')

% return

% identify VB volatile learner
dim = struct('n',2*5,'n_theta',3,'n_phi',2);
priors.muPhi = zeros(dim.n_phi,1);
priors.muTheta = [0;-2;0];
priors.muX0 = x0;
priors.SigmaPhi = 1e0*eye(dim.n_phi);
priors.SigmaTheta = 1e0*eye(dim.n_theta);
priors.SigmaX0 = 0e0*eye(dim.n);
priors.a_alpha = Inf;
priors.b_alpha = 0;
options.priors = priors;
[posterior,out] = VBA_NLStateSpaceModel(y,u,f_fname,g_fname,dim,options);

displayResults(posterior,out,y,x,x0,theta,phi,Inf,Inf)


[ha,hf] = unwrapVBvolatileOTO(posterior,out);
set(hf,'name','estimated volatile VB learner')




