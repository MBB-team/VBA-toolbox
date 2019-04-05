function [data, posterior, out] = demo_Kalman4bandits (data)
% // VBA toolbox //////////////////////////////////////////////////////////
%
% [data, posterior, out] = demo_Kalman4bandits (data)
%
% Demo of a multi-armed bandit task as described eg. in:
% Daw et al. (2006) Cortical substrates for exploratory decisions in
% humans. Nature 441(7095): 876-879
%
% If no inputs are given, the demo will generate a new experimental design 
% and simulate artificial behavioural responses.
%
% IN:
%   - data: optional structure with the following optional fields
%       - nTrials: number of trials {300}
%       - nBandits: number of bandits {4}
%       - payoff: nBandits X nTrials matrix specifying the outcome
%           associated to each bandit for each trial. If not specified, a new
%           payoff matrix will be generated using a random walk (see Daw et 
%           al. 2006 for details).
%       - choices: nBandits X nTrials matrix indicating for each trial and 
%           each bandit if the bandit has been chosen (1) or not (0). If not 
%           specified, artificial data will be simulated.
%
% OUT:
%   - data: same strucutre as the input with defaults filled in.
%   - posterior, out: results of the inversion 
%
% /////////////////////////////////////////////////////////////////////////

    % check inputs
    % =========================================================================
    if nargin == 0
        data = struct ();
    end

    % if no design specified, generate a new one
    % -------------------------------------------------------------------------
    if ~ isfield (data, 'payoff')
        data = VBA_check_struct (data, ...
            'nBandits', 4, ...
            'T', 300 ...
            );
        data = generateDesign (data);
    else
        [data.nBandits, data.T] = size (data.payoff);
    end

    % if no behavioural data given, simulate them
    % -------------------------------------------------------------------------
    if ~ isfield (data, 'choices')
        data = simulate4bandits (data);
    end

    % prepare data for model inversion
    % =========================================================================
    % observations
    y = data.choices;

    % inputs
    u = [ nan(data.nBandits, 1), data.choices ;  % previous choice
          nan, data.payoff(data.choices == 1)'] ; % previous feedback
    % note: the first input is a nan as there is no "previous trial" for the
    % first trial of the experiment

    % specify model
    % =========================================================================
    % observation and evolution functions
    % -------------------------------------------------------------------------
    f_fname = @f_Kalman4bandits; % evolution function (Kalman filter)
    g_fname = @g_Kalman4bandits; % observation function (softmax mapping)

    % model dimensions
    % -------------------------------------------------------------------------
    dim = struct( ...
        'n', 2 * data.nBandits, ... number of hidden states 
        'n_theta', 4, ... number of evolution parameters 
        'n_phi', 2 ... number of observation parameters 
       );

    % priors
    % -------------------------------------------------------------------------
    % initial state
    options.priors.muX0 = [50 * ones(data.nBandits, 1);  % muHat
                           5 * ones(data.nBandits ,1) ]; % sigmaHat

    options.priors.SigmaX0 = blkdiag( ...
        10 * eye (data.nBandits), ... % muHat
        1 * eye (data.nBandits)); % sigmaHat

    % evolution: lambda, theta, sigma_d, sigma_o
    options.priors.muTheta = [1; 50; 3; 4];
    options.priors.SigmaTheta = diag([1, 20, 1, 0]);

    % observation: beta, phi
    options.priors.muPhi = [log(0.1); 0];
    options.priors.SigmaPhi = diag([0.5, 0]);

    % task dimensions
    in = struct ('nBandits', data.nBandits);
    options.inF = in;
    options.inG = in;

    % options for the simulation
    % -------------------------------------------------------------------------
    % fitting categorical data
    options.sources.type = 2;

    % invert model
    % =========================================================================
    [posterior, out] = VBA_NLStateSpaceModel(y, u, f_fname, g_fname, dim, options);

end

% #########################################################################
% subroutines
% #########################################################################

% =========================================================================
% Task design generation. Draw random bandit/outcome pairings using a 
% gaussian random walk + noise, as described in Daw et al. 2006.
% =========================================================================
function data = generateDesign (data)
    % generative parameters
    % -------------------------------------------------------------------------
    % diffusion process parameters
    lambda = 0.9836;
    theta = 50;
    sigma_d = 2.8;
    % payoff variance
    sigma_o = 4;
        
    % loop over designs untill a satisfying one is found
    % -------------------------------------------------------------------------
    ok = false;
    while ~ ok
        % generate random walk
        mu(:, 1) = randi (100, data.nBandits, 1);
        for t = 2 : data.T
            % random step
            nu = VBA_random ('Gaussian', 0, sigma_d^2, data.nBandits, 1);
            % diffusion with decay
            mu(:, t) = lambda * mu(:, t - 1) + (1 - lambda) * theta + nu;    
        end

        % generate payoff by adding random noise
        payoff = round (mu + VBA_random ('Gaussian', 0, sigma_o^2, data.nBandits, data.T));

        % test if looks like a good design
        [~,mIdx] = max(payoff);
        counts = histcounts(mIdx,1:5);
        ok = VBA_isInRange(payoff, [0 100]) ... % outome between 0 and 100
             && all(counts > data.T/(2*data.nBandits)); % each bandit is the best at some point
    end
    
    % store parameters and payoff matrix
    % -------------------------------------------------------------------------
    data.design.parameters.lambda = lambda;
    data.design.parameters.theta = theta;
    data.design.parameters.sigma_d = sigma_d;
    data.design.parameters.sigma_o = sigma_o;
    data.design.mu = mu;
    data.payoff = payoff;
end

% =========================================================================
% Simulation routine. Simulate learning and decision of an agent
% implementing a Kalman filter belief update and a softmax decision rule,
% with parameters equal to the average subject.
% =========================================================================
function data = simulate4bandits (data)
    % create feedback structure for the simulation with VBA    
    % -------------------------------------------------------------------------
    % feeback rule: pick outcome from the payoff matrix
    h_feedback = @(yt,t,in) data.payoff(yt == 1, t); 

    % feedback structure
    fb = struct( ...
        'h_fname', h_feedback, ... % feedback function  
        'indy', 1 : data.nBandits, ... % where to store simulated choice
        'indfb', data.nBandits + 1, ... % where to store simulated feedback
        'inH', struct () ...
       );

    % define parameteters of the simulated agent    
    % -------------------------------------------------------------------------
    % cf Daw et al. 20016, Supplementary table 2

    % evolution parameters
    lambda = 0.924;
    theta = 50.5;
    sigma_d = 51.3;
    sigma_o = 4;

    pEvol = [lambda; theta; sigma_d; sigma_o];

    % observation parameters
    beta = 0.112;
    phi = 0;

    pObs = [log(beta); phi];

    % initial state
    x0 = [85.7 * ones(data.nBandits, 1);
          4.61 * ones(data.nBandits, 1) ];

    % pass on task dimension
    in = struct('nBandits', data.nBandits);
    options.inF = in;
    options.inG = in;

    % options for the simulation
    % -------------------------------------------------------------------------
    % number of trials
    n_t = data.T; 

    % fitting categorical data
    options.sources.type = 2;

    % simulate choices
    % -------------------------------------------------------------------------
    u = nan(data.nBandits+1, n_t);

    [y,x,~,~,~,u] = VBA_simulate ( ...
        n_t, ... number of trials
        @f_Kalman4bandits, ... evolution function
        @g_Kalman4bandits, ... observation function
        pEvol, ... evolution parameters (learning rate)
        pObs, ... observation parameters,
        u, ... dummy inputs
        Inf, Inf, ... deterministic evolution and observation
        options, ... options
        x0, ... initial state
        fb ... feedback rule
       );

    % store simulated choices, feedbacks, and parameters used
    % -------------------------------------------------------------------------
    data.choices = y;
    data.simulation.feedbacks = u(end, :);
    data.simulation.muHat = x(1 : data.nBandits, :);
    data.simulation.sigmaHat = x(data.nBandits + (1 : data.nBandits), :);
    data.simulation.parameters.lambda = lambda;
    data.simulation.parameters.theta = theta;
    data.simulation.parameters.sigma_d = sigma_d;
    data.simulation.parameters.sigma_o = sigma_o;
end
  

