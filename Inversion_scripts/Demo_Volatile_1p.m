%----------------------------------------------------------------------
% demo probability learning (from Mathys, Daunizeau et al. 2010
% Task : a sequence of binary variable has to be predicted
% Model returns his prediction at each trial before he learns from the actual real outcome.
%
% -- Reference :
% The model is based on the model used in :
% "A Bayesian foundation for individual learning under uncertainty"
% Christoph Mathys, Jean Daunizeau
%  doi:  10.3389/fnhum.2011.00039
%----------------------------------------------------------------------

close all
clear variables
clc

%----------------------------------------------------------------------
% DATA SIMULATION (for a single session)
%----------------------------------------------------------------------

% evolution, observation and feedback functions
f_fname = @f_VBvolatile_1p; % evolution function for the learning of probability
g_fname = @g_Mathys; % softmax decision based on probabilities
h_fname = @h_identity; % feedback is the actual correct choice 

% allocate feedback struture for simulations
u02 = [rand(1,50)<0.2]; 
fb.inH.u0 = repmat([u02,~u02],1,6); % definition of the binary time-series to be predicted
fb.h_fname = h_fname;
fb.indy = 2; % where to write subject choice in vector u
fb.indfb = 1; % where to write subject choice 

% simulation parameters % See  Mathys, Daunizeau et al. 2010 for a
% detailled description of parameters
theta = [1;-4;-1];
inF.lev2 = 1; % remove 3rd level (volatility learning)
inF.kaub = 1.4;
inF.thub = 1;
inF.rf = -1;
phi = []; 
inG.respmod = 'taylor';

% choose initial conditions
x0 = [0.5;0;0;1;log(4)];
u = zeros(2,size(fb.inH.u0,2)+1);

dim = struct('n',5,... % number of hidden states in the probability learning model
             'n_theta',3,...
             'n_phi',0);

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


% Simulate model and format data for further inversion
[y,x,x0,eta,e,u] = simulateNLSS_fb(length(u),f_fname,g_fname,theta,phi,u,Inf,Inf,options,x0,fb);


figure
hold on
plot(y-e,'r')
plot(y,'kx')
legend({'p(y=1|theta,phi,m)','binomial data samples'})
plot(fb.inH.u0,'go')
getSubplots
xlabel('time trials')

figure
hold on
plot(x')
title('hidden states')
xlabel('time trials')

%%


%----------------------------------------------------------------------
% MODEL INVERSION
%----------------------------------------------------------------------

% Model inversion
options.isYout = zeros(1,size(y,2));
[posterior,out] = VBA_NLStateSpaceModel(y,u,f_fname,g_fname,dim,options);
displayResults(posterior,out,y,x,x0,theta,phi,Inf,Inf)


