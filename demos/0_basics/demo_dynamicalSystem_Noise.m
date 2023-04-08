function [posterior, out] = demo_dynamicalSystem_Noise ()
% // VBA toolbox //////////////////////////////////////////////////////////
%
% [posterior, out] = demo_dynamicalSystem_Noise ()
% Demo of dynamical system inversion with two different
% methods of setting the observation noise prior.
%
% This demo provides a simple example of the estimation for a simple 
% dynamical system using informed priors on the observation noise.
%
% M. Eichenlaub 04/12/2019
% 
% /////////////////////////////////////////////////////////////////////////

%% Define the model
% =========================================================================

% Description of the dynamics
% -------------------------------------------------------------------------
% evolution is a parametric convolution of the input with one hidden state
f_fname = @f_alpha;
% - define size of a time step
dt = 1e-1; % in sec
% - store in the structure that will be passed to the evolution function
options.inF.dt = dt;

% observations are simply noisy mappings of the states (y = x + noise)
g_fname = @g_Id;
% however only the first state is observable
options.inG.ind = 1;

% Dimensions of the model
% -------------------------------------------------------------------------
% evolution parametesr
dim.n_theta = 2;
% observation parameters
dim.n_phi = 0;
% number of states
dim.n = 2;

%% Simulate data
% =========================================================================

% Define design, here inputs to be convolved
% -------------------------------------------------------------------------
% number of observations
n_t = 2e1 / dt; 
% baseline
u = zeros (1, n_t);
% random pulses
nPulses = 4;
u(randperm (n_t, nPulses)) = 1;

% Parameters of the model to be simulated
% -------------------------------------------------------------------------
% evolution parameters
theta = [1.1 / dt; 0.1];
% observation parameters
phi = [];
% initial state
x0 = [0; 0];

% state noise precision
alpha = Inf; % deterministic  

% observation precision
% Setting true mean the observation noise standard deviation (SD).
m_true = 0.01; 

sigma = 1/(m_true)^2; % Caluclating assocoited precison

% Simulate
% -------------------------------------------------------------------------
[y, x] = VBA_simulate (n_t, f_fname, g_fname, theta, phi, u, alpha, sigma, options, x0);

% Display
% -------------------------------------------------------------------------
plotSimulation (u, x, y, n_t, dt);

%% Perform model estimation using VBA_guessHyperpriors
% =========================================================================

% Option 1: 
% Call VBA_guessHyperpriors(y) which estimates the observation noise from
% the data directly.

[a,b] = VBA_guessHyperpriors(y);
options.priors.a_sigma = a;
options.priors.b_sigma = b;

% call VBA invesrion routine
% -------------------------------------------------------------------------
options.verbose = false;
[posterior,out] = VBA_NLStateSpaceModel(y,u,f_fname,g_fname,dim,options);

% Get Measurement noise results

% Caluculate mean and SD from observation noise SD prior
[m_pr,s_pr] = VBA_Convert_ab(options.priors.a_sigma,options.priors.b_sigma);

% Caluculate mean and SD from observation noise SD posterior results
[m,s] = VBA_Convert_ab(posterior.a_sigma,posterior.b_sigma);

fprintf('\nResults using VBA_guessHyperpriors:\n')
fprintf('Prior Noise SD: %.3f +/- %.3f\n',m_pr,s_pr);
fprintf('Posterioir Noise SD: %.3f +/- %.3f\n',m,s);
fprintf('True Noise SD: %.3f\n',m_true);


%% Perform model estimation using VBA_Create_NoisePrior
% =========================================================================

% Option 2: 
% Call VBA_Create_NoisePrior(m,s) which allows the specification of the
% prior distribution of the observation noise SD directly, using a mean 
% and SD. This is useful when prior knowledge on the observation error is
% available.

[a,b] = VBA_Create_NoisePrior(0.1,0.1); 
options.priors.a_sigma = a;
options.priors.b_sigma = b;

% call VBA invesrion routine
% -------------------------------------------------------------------------
options.verbose = false;
[posterior,out] = VBA_NLStateSpaceModel(y,u,f_fname,g_fname,dim,options);

% Get observation noise results

% Caluculate mean and SD from observation noise SD prior. This recovers
% the values set before for the prior
[m,s] = VBA_Convert_ab(posterior.a_sigma,posterior.b_sigma);

% Caluculate mean and SD from observation noise SD posterior results
[m_pr,s_pr] = VBA_Convert_ab(options.priors.a_sigma,options.priors.b_sigma);

fprintf('\nResults using VBA_Create_NoisePrior:\n')
fprintf('Prior Noise SD: %.3f +/- %.3f\n',m_pr,s_pr);
fprintf('Posterioir Noise SD: %.3f +/- %.3f\n',m,s);
fprintf('True Noise SD: %.3f\n',m_true);

end

%% ########################################################################
% Display functions
% ########################################################################
function plotSimulation (u, x, y, n_t, dt)
    hf = figure ('color', 'w', 'Name', 'demo_dynamicalSystem: simulation');

    % x
    timeline = (0 : n_t - 1) * dt;

    % plot inputs
    ha(1) = subplot (3, 1, 1, 'parent', hf);
    plot (ha(1), timeline, u);
    xlabel (ha(1), 'time (sec)');
    ylabel (ha(1), 'u(t)');
    title (ha(1), 'input');

    % plot state trajectory
    ha(2) = subplot (3, 1, 2, 'parent', hf);
    plot (ha(2), timeline, x);
    xlabel (ha(2), 'time (sec)');
    ylabel (ha(2), 'x(t)');
    title (ha(2), 'state');

    % plot observations
    ha(3) = subplot (3, 1, 3, 'parent', hf);
    plot (ha(3), timeline, y, '.');
    xlabel (ha(3), 'time (sec)');
    ylabel (ha(3), 'y(t)');
    title (ha(3), 'observation');

    % make pretty
    set (ha, 'box', 'off', 'ygrid', 'on', 'ylim', [- 0.2, 1.2]);
end

