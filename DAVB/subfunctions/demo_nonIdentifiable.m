% Demo for non-identifiable parameters

close all
clear variables

%---- simulate noisy gbf ----%

% Choose basic settings for simulations
sigma = 1e1;            % precision 
phi = [5;-1];       % observation parameters
g_fname = @g_NI;        % observation function
inG.n = 1e2;             % sample size
inG.X = ones(inG.n,1);
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
priors.muPhi = zeros(2,1);         % prior mean on observation params
priors.SigmaPhi = 1e1*eye(2); % prior covariance on observation params
% priors.SigmaPhi(2,2) = 0;
priors.a_sigma = 1e0;             % Jeffrey's prior
priors.b_sigma = 1e0;             % Jeffrey's prior
% Build options structure
% options.checkGrads = 1;
options.priors = priors;        % include priors in options structure
options.inG = inG;              % input structure (grid)
options.GnFigs = 0;             % disable annoying figures
options.verbose = 1;
dim.n_phi = 2;                  % nb of observation parameters
dim.n_theta = 0;                % nb of evolution parameters
dim.n=0;                        % nb of hidden states
% Call inversion routine
[posterior,out] = VBA_NLStateSpaceModel(y,[],[],g_fname,dim,options);


%---- Display results ----%
displayResults(posterior,out,y,[],[],[],phi,[],sigma)

