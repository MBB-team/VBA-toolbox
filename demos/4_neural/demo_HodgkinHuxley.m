function [posterior, out] = demo_HodgkinHuxley ()
% // VBA toolbox //////////////////////////////////////////////////////////
%
% [posterior, out] = demo_HodgkinHuxley ()
% Demo of Q-learning simulation and inference
%
% Demo of the Hodgin-Huxley neuronal model.
%
% Background:
% ~~~~~~~~~~~

% Input current (u) are sent to a neuron that responds according to the
% Hodgkin-Huxley model. In brief, an action potential is generated when the
% membrane depolarization reaches a critical threshold. An AP
% approximately corresponds to a 80mV depolarization. 
% This demo simulates the response of such a neuron, and then inverts the
% model. Emphasis is put on the identifiability of model parameters 
% (e.g. K/Na conductances).
%
% /////////////////////////////////////////////////////////////////////////

% Basic settings for simulations
% =========================================================================
% Sampling rate 
delta_t = 1 / 10; % 10Hz
% Recording duration
n_t = 500; 

% number of integration timesteps between two observation
options.decim = 3; % smaller means faster but high risk of numerical errors

% specify model
% =========================================================================
% evolution function
% -------------------------------------------------------------------------
f_fname = @f_HH;

% stepsize of the numerical integration
options.inF.delta_t = delta_t / options.decim;

% observation function
% -------------------------------------------------------------------------
g_fname = @g_Id;

% only membrane potential is measured
options.inG.ind = 1;

% Simulation
% =========================================================================
% input currents from random spikes
pSpike = 0.05;
spikeAmp = 100;
u = spikeAmp * VBA_random ('Bernoulli', pSpike, 1, n_t);

% Parameters of the simulation
steadyState = [0; - 2.8; - 0.75; 0.4]; % asymptotic state with 0 input
theta = 0.1 * randn (4, 1); % evolution parameters, 0 = default HH params
phi = [];
alpha = Inf;
sigma = 1e2;

% simulate
[y, x, x0, eta, e] = VBA_simulate (n_t,f_fname,g_fname,theta,phi,u,alpha,sigma,options,steadyState);

% display
displayHH (u, x);

% Inversion: parameter recovery
% =========================================================================
% dimensions of the problem
dim.n = 4;
dim.n_theta = 4;
dim.n_phi = 0;

% priors
options.priors.muX0 = steadyState;

% estimation routine
[posterior,out] = VBA_NLStateSpaceModel(y,u,f_fname,g_fname,dim,options);

% display
displayResults(posterior,out,y,x,x0,theta,phi,alpha,sigma);

end

%% subfunctions
% #########################################################################
function displayHH (u, x)

hf = figure('color',[1 1 1]);
ha = subplot(3,1,1,'parent',hf,'nextplot','add');
plot(ha,u')
title(ha,'input current')
xlabel('time')
ha = subplot(3,1,2,'parent',hf,'nextplot','add');
col = getColors(2);
m = VBA_sigmoid(x(2,:));
n = VBA_sigmoid(x(3,:));
h = VBA_sigmoid(x(4,:));
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
end


