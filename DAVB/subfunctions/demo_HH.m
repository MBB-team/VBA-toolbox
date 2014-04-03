% Demo for Hodgkin-Huxley PA propagation
% Short bursts of depolarizing input current (u) are sent to a neuron that
% responds according to Hodgkin-Huxley model. In brief, an AP is generated
% if the membrane depolarization reaches a critical threshold. NB: an AP
% approximately corresponds to a 80mV depolarization. This script simulates
% the response of such a neuron, and then inverts the model. Emphasis is
% put on the identifiability of model parameters (e.g. K/Na conductances).

clear variables
close all

% Choose basic settings for simulations
n_t = 1e3;
delta_t = 1e-1;         % 10Hz sampling rate
f_fname = @f_HH;
g_fname = @g_Id;

u       = 4e1*(randn(1,n_t)>1.1);
% figure,plot(u)

% Build options structure for temporal integration of SDE
inF.delta_t = delta_t;
inG.ind         = 1;
options.inF     = inF;
options.inG     = inG;
options.decim = 1;


% Parameters of the simulation
alpha   = Inf;
sigma   = 1e2;
phi     = [];
x0 = [0;-3;-3/4;1/3];
theta = [1;0*ones(3,1)];

% Build time series of hidden states and observations
[y,x,x0,eta,e] = simulateNLSS(n_t,f_fname,g_fname,theta,phi,u,alpha,sigma,options,x0);

% display time series of hidden states and observations
displaySimulations(y,x,eta,e)

hf = figure('color',[1 1 1]);
ha = subplot(3,1,1,'parent',hf,'nextplot','add');
plot(ha,u')
title(ha,'input current')
xlabel('time')
ha = subplot(3,1,2,'parent',hf,'nextplot','add');
col = getColors(2);
m = sigm(x(2,:));
n = sigm(x(3,:));
h = sigm(x(4,:));
p(1,:) = m.^3.*h;
p(2,:) = n.^4;
for i=1:2
    plot(ha,p(i,:),'color',col(i,:));
end
title(ha,'ion channels opening probabilities')
xlabel('time')
legend({'Na','K'})
ha = subplot(3,1,3,'parent',hf,'nextplot','add');
plot(ha,x(1,:))
title(ha,'output membrane depolarization (mV)')
xlabel('time')

% return

dim.n = 4;
dim.n_phi = 0;
dim.n_theta = 4;

priors.iQx = cell(n_t,1);
for t=1:n_t
    priors.iQx{t} = 1e2*eye(4);
    priors.iQx{t}(1,1) = 1e0;
end
priors.a_alpha = Inf;
priors.b_alpha = 0;
options.priors = priors;

[posterior,out] = VBA_NLStateSpaceModel(y,u,f_fname,g_fname,dim,options);


% Display results
displayResults(posterior,out,y,x,x0,theta,phi,alpha,sigma)


