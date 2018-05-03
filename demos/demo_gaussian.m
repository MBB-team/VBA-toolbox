% Basic demo for fitting mixture of Gaussians
% This demo first simulates a noisy MoG density profile, whose mode
% position and width are then estimated using VBA.

close all
clear all
clc

%---- simulate noisy data
% Choose basic settings for simulations
sigma = 1e-1;            % precision 
phi = [5;-1;1;3];         % observation parameters
g_fname = @g_Gaussian; % observation function
inG.grid = -10:0.1:10;  % grid on which the gbf is evaluated
inG.input = zeros(length(inG.grid),1);
inG.input(50) = 1;
% Build simulated observations 
[gx] = feval(g_fname,[],phi,[],inG);
y = gx + sqrt(sigma.^-1)*randn(size(gx));
% display time series of hidden states and observations
figure,
plot(y','ro')
hold on;
plot(gx')


%---- Invert model on simulated data
% Build priors structure
priors.muPhi = zeros(4,1); % prior mean on observation params
priors.SigmaPhi = 1e0*eye(4); % prior covariance on observation params
priors.a_sigma = 1; % Jeffrey's prior
priors.b_sigma = 1; % Jeffrey's prior
options.priors = priors; % include priors in options structure
options.inG = inG; % input structure (grid)
dim.n_phi = 4; % nb of observation parameters
dim.n_theta = 0; % nb of evolution parameters
dim.n = 0;  % nb of hidden states
% Call inversion routine
[posterior,out] = VBA_NLStateSpaceModel(y,[],[],g_fname,dim,options);


%---- Display results
displayResults(posterior,out,y,[],[],[],phi,[],sigma);

