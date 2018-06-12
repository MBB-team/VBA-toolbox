function [posterior, out]=demo_Qlearning(choices, feedbacks)
% // VBA toolbox //////////////////////////////////////////////////////////
%
% [posterior, out] = demo_Qlearning([choices, feedbacks])
% Demo of Q-learning simulation and inference
%
% This is a simple example of reinforcement learning algorithm.
%
% Background:
% ~~~~~~~~~~~
% In psychological terms, motivation can be defined as the set of processes
% that generate goals and thus determine behaviour. A goal is nothing else
% than a 'state of affairs', to which people attribute (subjective) value.
% Empirically speaking, one can access these values by many means,
% including subjective verbal report or decision making. These measures
% have been used to demonstrate how value change as people learn a new
% operant response. This is predicted by reinforcement learning theories,
% which essentially relate behavioural response frequency to reward. In
% this context, value is expected reward, and it changes in proportion to
% the agent's prediction error, i.e. the difference between actual and
% expected reward.
%
% /////////////////////////////////////////////////////////////////////////

% check inputs
% =========================================================================

switch nargin
    case 0
        fprintf('No inputs provided, generating simulated behavior...\n\n');
        [choices, feedbacks, simulation]=demo_QlearningSimulation();
    case 2
        fprintf('Performing inversion of provided behaviour...\n\n');
    otherwise
        error('*** Wrong number of arguments.')
end
        
% reformat data
% =========================================================================
% observations
y = choices;
% inputs
u = [ nan, choices(1:end-1)   ;  % previous choice
      nan, feedbacks(1:end-1) ]; % previous feedback

% specify model
% =========================================================================
f_fname = @f_Qlearning; % evolution function (Q-learning)
g_fname = @g_QLearning; % observation function (softmax mapping)

% provide dimensions
dim = struct( ...
    'n', 2, ... number of hidden states (2 Q-values)
    'n_theta', 1, ... number of evolution parameters (1: learning rate)
    'n_phi', 1 ... number of observation parameters (1: temperature)
   );
    
% options for the simulation
% -------------------------------------------------------------------------
% use the default priors except for the initial state
options.priors.muX0 = [0.5; 0.5];
options.priors.SigmaX0 = 0.1 * eye(2);

% options for the simulation
% -------------------------------------------------------------------------
% number of trials
n_t = numel(choices); 
% fitting binary data
options.sources.type = 1;
% Normally, the expected first observation is g(x1), ie. after
% a first iteratition x1 = f(x0, u0). The skipf flag will prevent this evolution
% and thus set x1 = x0
options.skipf = [1 zeros(1,n_t-1)];

% invert model
% =========================================================================
[posterior, out] = VBA_NLStateSpaceModel(y, u, f_fname, g_fname, dim, options);

% display estimated parameters:
% -------------------------------------------------------------------------
fprintf('=============================================================\n');
fprintf('\nEstimated parameters: \n');
fprintf('  - learning rate: %3.2f\n', VBA_sigmoid(posterior.muTheta));
fprintf('  - inverse temp.: %3.2f\n\n', exp(posterior.muPhi));
fprintf('=============================================================\n');

% invert model
% =========================================================================
if exist('simulation','var') % used simulated data from demo_QlearningSimulation
    displayResults( ...
        posterior, ...
        out, ...
        choices, ...
        simulation.state, ...
        simulation.initial, ...
        simulation.evolution, ...
        simulation.observation, ...
        Inf, Inf ...
     );
end

