% Demo for model evidence penalty of parameter redundancy

close all
clear variables

%---- simulate noisy gbf ----%

% Choose basic settings for simulations
sigma = 1e1;            % precision 
phi = [1];         % observation parameters
g_fname = @g_dummy; % observation function
inG.X = randn(64,1);  % grid on which the gbf is evaluated
% Build simulated observations 
[gx] = feval(g_fname,[],phi,[],inG);
y = gx + sqrt(sigma.^-1)*randn(size(gx));
% display time series of hidden states and observations
figure,
plot(y','ro')
hold on;
plot(gx')
% disp('--paused--')
% pause


%---- Invert gbf on simulated data ----%

% Build priors structure
priors.muPhi = [0];         % prior mean on observation params
priors.SigmaPhi = 1e4*eye(1); % prior covariance on observation params            % Jeffrey's prior
% Build options structure
options.priors = priors;        % include priors in options structure
options.inG = inG;              % input structure (grid)
options.GnFigs = 0;             % disable annoying figures
options.verbose = 1;
dim.n_phi = 1;                  % nb of observation parameters
dim.n_theta = 0;                % nb of evolution parameters
dim.n=0;                        % nb of hidden states
% Call inversion routine
[p1,o1] = VBA_NLStateSpaceModel(y,[],[],g_fname,dim,options);

priors.muPhi = [0;0];         % prior mean on observation params
priors.SigmaPhi = 1e4*eye(2);
options.priors = priors;  
dim.n_phi = 2; 
[p2,o2] = VBA_NLStateSpaceModel(y,[],[],g_fname,dim,options);

o1.F - o2.F

