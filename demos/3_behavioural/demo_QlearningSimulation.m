function [choices, feedbacks, simulation]=demo_QlearningSimulation()
% // VBA toolbox //////////////////////////////////////////////////////////
%
% [choices, feedbacks] = demo_QlearningSimulation()
% Demo of Q-learning simulation. This demonstrate in particular how to use
% the toolbox to simulate data with a dynamic feedack (behaviour at times t
% can serve directly or indirectly as an input at time t+1)
%
% This is based on simple example of reinforcement learning algorithm.
% This demo first simulates 150 choices of a Q-learning agent when faced to
% a sequence binary alternatives: 
% - at each trial, the agent choose one of two actions
% - the action is rewarded according to the probabilistic contingency rule
% - every 25 trials, the contingencies are reversed.
%
% /////////////////////////////////////////////////////////////////////////

f_fname = @f_Qlearning; % evolution function (Q-learning)
g_fname = @g_QLearning; % observation function (softmax mapping)

% Create the feedback rule for the simulation
% =========================================================================

% define which action should be rewarded at each trial (contingencies)
% -------------------------------------------------------------------------
% probability of a positive reward following a 'correct' action 
    probRewardGood = 75/100;
% draw 25 random feedbacks
contBloc = +(rand(1,25) < probRewardGood); 
% create 6 blocs with reversals
contingencies = [contBloc, 1-contBloc, ...
                 contBloc, 1-contBloc, ...
                 contBloc, 1-contBloc] ;
        
% create feedback structure for the simulation with VBA    
% -------------------------------------------------------------------------
% feedback function. Return 1 if action follow contingencies.
h_feedback = @(yt,t,in) +(yt == contingencies(t));
% feedback structure for the VBA
fb = struct( ...
    'h_fname', h_feedback, ... % feedback function  
    'indy', 1, ... % where to store simulated choice
    'indfb', 2, ... % where to store simulated feedback
    'inH', struct() ...
   );

% Simulate choices for the given feedback rule
% =========================================================================

% define parameteters of the simulated agent    
% -------------------------------------------------------------------------
% learning rate
theta = VBA_sigmoid(0.65,'inverse',true); % 0.65, once sigm transformed
% inverse temperature 
phi = log(2.5); % will be exp transformed
% initial state
x0 = [.5; .5];

% options for the simulation
% -------------------------------------------------------------------------
% number of trials
n_t = numel(contingencies); 
% fitting binary data
options.sources.type = 1;
% Normally, the expected first observation (choice) is g(x1), ie. after
% a first iteratition x1 = f(x0). The skipf flag will prevent this evolution
% and thus set x1 = x0
options.skipf = [1 zeros(1,n_t)];

% simulate choices
% -------------------------------------------------------------------------
[y,x,x0,eta,e,u] = VBA_simulate ( ...
    n_t+1, ... number of trials
    f_fname, ... evolution function
    g_fname, ... observation function
    theta, ... evolution parameters (learning rate)
    phi, ... observation parameters,
    nan(2,n_t), ... dummy inputs
    Inf, Inf, ... deterministic evolution and observation
    options, ... options
    x0, ... initial state
    fb ... feedback rule
   );

% plot simulated choices
% -------------------------------------------------------------------------
hf = figure( ...
    'name', 'Simulated Q-learning behaviour', ...
    'color','w' ...
   );
ha = axes('parent',hf,'nextplot','add');
plot(ha,y,'kx')
plot(ha,y-e,'r')
legend(ha,{'y: agent''s choices','p(y=1|theta,phi,m): behavioural tendency'})

% Return simulated choices, feedbacks, and parameters used for the
% simulation
% =========================================================================
choices = u(1,2:end);
feedbacks = u(2,2:end);
simulation = struct( ...
    'state', x(:,1:end-1), ...
    'initial', x0, ...
    'evolution', theta, ...
    'observation', phi ...
    );


