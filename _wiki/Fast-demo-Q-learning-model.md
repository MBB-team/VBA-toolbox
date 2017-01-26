---
title: "A fast demo: Q-learning"
---
* Will be replaced with the ToC, excluding the "Contents" header
{:toc}

Below, we copied and pasted the content and graphical output of a demonstration script (`demo_Qlearning.m`) that replicates a typical learning study in behavioural economics and/or experimental psychology.
The task is the so-called two-armed bandit problem, which captures the essence of operant learning. The demo uses the so-called [Q-learning model]({{ site.baseurl }}/wiki/Structure-of-VBA's-generative-model) which predicts how people change their behavioural response according to the consequences of their actions (e.g., reward/punishment feedback).



```matlab
% Q-learning demo
% In psychological terms, motivation can be defined as the set of processes
% that generate goals and thus determine behaviour. A goal is nothing else
% than a “state of affairs”, to which people attribute (subjective) value.
% Empirically speaking, one can access these values by many means,
% including subjective verbal report or decision making. These measures
% have been used to demonstrate how value change as people learn a new
% operant response. This is predicted by reinforcement learning theories,
% which essentially relate behavioural response frequency to reward. In
% this context, value is expected reward, and it changes in proportion to
% the agent prediction error, i.e. the difference between actual and
% expected reward.
% This demo simulates a number of sequences of choices of a Q-learning
% agent, which is a simple example of reinforcement learning algiorithms.
% We then invert the model using VBA. Finally, we perform a Volterra
% decomposition of hidden states dynamics onto a set of appropriately
% chosen basis functions (here: the agent chosen action, and the winning
% action). This diagnostic analysis allows one to identify the hidden
% states impulse response to experimentally controlled inputs to the
% system.

close all
clear variables
clc


f_fname = @f_Qlearn2; % evolution function (Q-learning)
g_fname = @g_softmax; % observation function (softmax mapping)
h_fname = @h_randOutcome; % feedback function (reward schedule)

% allocate feedback struture for simulations
fb.inH.er = 1;
fb.inH.vr = 0;
fb.h_fname = h_fname;
fb.indy = 1;
fb.indfb = 2;
u0 = [randn(1,25)>-0.25]; % possible feedbacks
fb.inH.u0 = [u0,~u0,u0,~u0,u0,~u0]; % with reversals

% simulation parameters
theta = sigm(0.75,struct('INV',1)); % learning rate = 0.75
phi = log(2); % inverse temperature = 2
x0 = zeros(2,1);
n_t = size(fb.inH.u0,2)+1; % number of trials
options.binomial = 1;
options.verbose = 0;
options.skipf = zeros(1,n_t);
options.skipf(1) = 1; % apply identity mapping from x0 to x1.

% simulate Q-learner and plot choices
[y,x,x0,eta,e,u] = simulateNLSS_fb( ...
    n_t,f_fname,g_fname,theta,phi,zeros(2,n_t),Inf,Inf,options,x0,fb);
hf = figure('color',[1 1 1]);
ha = axes('parent',hf,'nextplot','add');
plot(ha,y,'kx')
plot(ha,y-e,'r')
legend(ha,{'y: agent''s choices','p(y=1|theta,phi,m): behavioural tendency'})


% VBA model inversion (given simulated choices)
dim = struct('n',2,'n_theta',1,'n_phi',1);
priors.a_alpha = Inf;
priors.b_alpha = 0;
options.priors = priors;
[posterior,out] = VBA_NLStateSpaceModel(y,u,f_fname,g_fname,dim,options);

% compare simulated and estimated model variables
displayResults(posterior,out,y,x,x0,theta,phi,Inf,Inf);


% perform Volterra decomposition
u1 = u(1,:); % own action
u3 = u(2,:); % feedback
u2 = zeros(size(u1)); % opponent action
u2(u3>0) = u1(u3>0);
u2(u3<0) = 1-u1(u3<0);
uu = 2*[u1;u2]-1;
o = out;
o.u = uu;
[kernels] = VBA_VolterraKernels(posterior,o,16);
o.diagnostics.kernels = kernels;
VBA_ReDisplay(posterior,o,1)

getSubplots
```

Below are graphical outputs of the demonstration script. First, let us focus on the simulated behavioural tendency and agent's choices:

![]({{ site.baseurl }}/images/wiki/demo1/demo1_1.jpg)

On can see how the agent's behavioural response changes according to the feedback he receives. In brief, the Q-learner is tracking the winning option (which effectively varies over time).

All VBA's graphical outputs (parameter estimation, model accuracy, inversion diagnostics, convergence, etc...), for this particular example, are described [here]({{ site.baseurl }}/wiki/VBA-graphical-output). Now let us check how accurate the model inversion was:

![]({{ site.baseurl }}/images/wiki/demo1/demo1_2.jpg)

> **Upper-left panel**: posterior mean (grey bar) and standard deviation (red errorbar) of the model's evolution parameter (theta = learning rate, here). The simulated parameter is shown using a green dot. **Upper-right panel**: same for the observation parameter (phi = behavioural log-temperature). **Middle-up-right panel**: same for initial conditions (options values). **Middle-down-left panel**: same for hidden states (2 option values plotted as a function of time or trials). **Middle-down-right panel**: same for observed data. **Lower-left panel**: estimated hidden states (y-axis) are plotted against simulated hidden states (x-axis). **Lower-left panel**: fitted data (y-axis) are plotted against simulated data (x-axis).

Recall the unknown parameters of the [Q-learning model]({{ site.baseurl }}/wiki/Structure-of-VBA's-generative-model): 1 evolution parameter (learning rate), 1 observation parameter (temperature) and 2 initial conditions (initial action values).
One can see that the posterior credible intervals (red errorbars) contain the simulated parameter values (green dots). In turn, estimated and simulated action values (as well as choices) dynamics tightly correlate with each other.

Of course, analysis of experimental data does not allow one to evaluate the accuracy of parameter estimation, as was done here using simulated data. However, VBA provides a number of [inversion diagnostics]({{ site.baseurl }}/wiki/VBA-output-structure) that are useful to review, such as the structure of model residuals, the parameter posterior correlation matrix, Volterra kernels, etc...

Note also that VBA also includes routines for performing other types of analyses. These include, but are not limited to:

- **clustering analyses**: in this type of analysis, data are assumed to be sampled from a mixture of Gaussian or Binomial distributions.
- **group-level Bayesian model comparisons**: In this type of analysis, models are treated as random effects that could differ between subjects and have a fixed (unknown) distribution in the population.

Finally, note that VBA already includes a large library of models for behavioural and neurophysiological data...

