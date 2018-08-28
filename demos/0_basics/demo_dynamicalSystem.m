
function [posterior, out] = demo_dynamicalSystem ()
% // VBA toolbox //////////////////////////////////////////////////////////
%
% [posterior, out] = demo_dynamicalSystem ()
% Demo of dynamical system simulation and inversion
%
% This demo provides a simple example of the simulation and the estimation
% for a simple dynamical system.
% As a comparison, this demo also implement a grid search that can be
% compared to the variational result.
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
n_t = 1e2 / dt; 
% baseline
u = zeros (1, n_t);
% random pulses
nPulses = 16;
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
sigma = 1e2; % no noise

% Simulate
% -------------------------------------------------------------------------
[y, x] = VBA_simulate (n_t, f_fname, g_fname, theta, phi, u, alpha, sigma, options, x0);

% Display
% -------------------------------------------------------------------------
plotSimulation (u, x, y, n_t, dt);

%% Perform model estimation, the VBA way
% =========================================================================
fprintf('Variational estimation... ');
tic;

% call VBA invesrion routine
% -------------------------------------------------------------------------
options.verbose = false;
[posterior,out] = VBA_NLStateSpaceModel(y,u,f_fname,g_fname,dim,options);

% find best parameter
% -------------------------------------------------------------------------
theta_VB = posterior.muTheta;

% chrono
% -------------------------------------------------------------------------
fprintf('done! (took %f seconds)\n', toc);

%% Perform model estimation, the old way
% =========================================================================
% for the sake of the example, we'll now try to recover the parameters from
% the data using a grid search over parameters and look for the set that 
% maximimizes the likelihood
fprintf('Performing grid search... ');
tic;

% define grid for parameter search
% -------------------------------------------------------------------------
N = 9;
g_theta1 = 1 / dt + linspace (- 4 , 4,  N);
g_theta2 = linspace (- 1, 1, N);

% Compute predictions and likelihoodfor each parameter set
% -------------------------------------------------------------------------
options.verbose = false;
for i1 = 1 : numel (g_theta1)
    for i2 = 1 : numel (g_theta2)
        theta_grid = [g_theta1(i1); g_theta2(i2)];
        sigma_test = Inf;
        g{i1, i2} = VBA_simulate (n_t, f_fname, g_fname, theta_grid, phi, u, alpha,sigma_test, options, x0);
        LL(i1, i2) = - sum ( (y - g{i1, i2}) .^ 2);
    end
end

% find best parameter
% -------------------------------------------------------------------------
[~, i1_ML, i2_ML] = VBA_maxMat(LL);
theta_ML = [g_theta1(i1_ML); g_theta2(i2_ML)];

% Display results
% -------------------------------------------------------------------------
plotGridSearch (y, g, g_theta1, g_theta2, LL, n_t, dt);

% chrono
% -------------------------------------------------------------------------
fprintf('done! (took %f seconds)\n', toc);

%% Show results
% =========================================================================
fprintf('Results of the estimation:\n');
disp (table (...
    [theta(1); theta_ML(1); theta_VB(1)], ...
    [theta(2); theta_ML(2); theta_VB(2)], ...
    'RowNames',{'trueValue','ML','VB'}, ...
    'VariableNames',{'theta_1','theta_2'}));

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

function plotGridSearch (y, g, g_theta1, g_theta2, LL, n_t, dt)
    hf = figure ('color', 'w', 'Name', 'demo_DynamicalSystem: grid search');

    % x
    timeline = (0 : n_t - 1) * dt;

    % Likelihood
    eLL = exp(LL/100);
    hb(1) = subplot (1, 2, 2, 'parent', hf);
    imagesc(eLL, 'Parent', hb(1));
    colorLims = [-0.15 1.1*VBA_maxMat(eLL)];
    set(hb(1),'XTick',1 : numel(g_theta2), 'XTickLabel', g_theta2);
    set(hb(1),'YTick',1 : numel(g_theta1), 'YTickLabel', g_theta1);
    set(hb(1),'Clim',colorLims);
    colorbar
    title (hb(1), 'Likelihood');
    xlabel (hb(1), 'theta_2')
    ylabel (hb(1), 'theta_1')

    % predictions about theta 
    hb(2) = subplot (1, 2, 1, 'parent', hf);
    xlabel (hb(2), 'time (sec)');
    ylabel (hb(2), 'observations / predictions');
    title (hb(2), 'trajectories');
    hold on;

    % make pretty
    set (hb(2), 'box', 'off', 'ygrid', 'on');

    % color predictions as function of likelihood
    xLL = linspace(colorLims(1),colorLims(2), 64);
    cLL = colormap(flipud(colormap('hot')));
    for i1 = 1 : numel (g_theta1)
        for i2 = 1 : numel (g_theta2)
        l = plot(hb(2),timeline,g{i1,i2});
        %l(2) = plot(hb(3),timeline,g{i1,i2}(2,:));
        set(l,'Color',interp1(xLL, cLL, eLL(i1,i2)));
        end
    end

    % add data
    plot (hb(2), timeline, y, '.', 'Color', [.6 .6 .6], 'MarkerSize', 6);
end
