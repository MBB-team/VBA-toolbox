---
title: "A fast demo: Q-learning"
---
* Will be replaced with the ToC, excluding the "Contents" header
{:toc}

Below, we copied and pasted the content and graphical output of a demonstration script (`demo_Qlearning.m`) that replicates a typical learning study in behavioural economics and/or experimental psychology.
The task is the so-called two-armed bandit problem, which captures the essence of operant learning. The demo uses the so-called [Q-learning model]({{ site.baseurl }}/wiki/Structure-of-VBA's-generative-model) which predicts how people change their behavioural response according to the consequences of their actions (e.g., reward/punishment feedback).


```matlab
f_fname = @f_Qlearn2; % evolution function (Q-learning)
g_fname = @g_softmax; % observation function (softmax mapping)
h_fname = @h_randOutcome; % feedback function (reward schedule)

% allocate feedback struture for simulations (see simulateNLSS_fb.m)
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
x0 = zeros(2,1); % initial action values = 0
n_t = size(fb.inH.u0,2)+1; % number of trials
options.binomial = 1; % choice data = 0 or 1
options.skipf = [1,zeros(1,n_t-1)]; % apply identity mapping from x0 to x1.

% simulate Q-learner and plot choices
[y,x,x0,eta,e,u] = simulateNLSS_fb(n_t,f_fname,g_fname,theta,phi,zeros(2,n_t),Inf,Inf,options,x0,fb);
hf = figure('color',[1 1 1]);
ha = axes('parent',hf,'nextplot','add');
plot(ha,y,'kx')
plot(ha,y-e,'r')
legend(ha,{'y: agent''s choices','p(y=1|theta,phi,m): behavioural tendency'})


% VBA model inversion (given simulated choices)
dim = struct('n',2,'n_theta',1,'n_phi',1); % 2 hidden states, 1 evolution param, 1 observation param
[posterior,out] = VBA_NLStateSpaceModel(y,u,f_fname,g_fname,dim,options);

% compare simulated and estimated model variables
displayResults(posterior,out,y,x,x0,theta,phi,Inf,Inf);
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

