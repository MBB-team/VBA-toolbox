% Demo for fitting mixture of Gaussians

close all
clear variables

%---- simulate noisy gbf ----%

% Choose basic settings for simulations
sigma = 1e1;            % precision 
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
% disp('--paused--')
% pause


%---- Invert gbf on simulated data ----%

% Build priors structure
priors.muPhi = zeros(4,1);         % prior mean on observation params
priors.SigmaPhi = 1e4*eye(4); % prior covariance on observation params
priors.a_sigma = 1;             % Jeffrey's prior
priors.b_sigma = 1;             % Jeffrey's prior
% priors.iQy{1} = 0*eye(length(gx));
% priors.iQy{1}(150,150) = 1;
% priors.iQy{1}(250,250) = 1;
% priors.iQy{1}(100,100) = 1;
% priors.iQy{1}(200,200) = 1;
% priors.iQy{1}(300,300) = 1;
% Build options structure
% options.checkGrads = 1;
options.Laplace = 1;
options.priors = priors;        % include priors in options structure
options.inG = inG;              % input structure (grid)
options.GnFigs = 0;             % disable annoying figures
options.verbose = 1;
dim.n_phi = 4;                  % nb of observation parameters
dim.n_theta = 0;                % nb of evolution parameters
dim.n=0;                        % nb of hidden states
% Call inversion routine
[posterior,out] = VBA_NLStateSpaceModel(y,[],[],g_fname,dim,options);


%---- Display results ----%
displayResults(posterior,out,y,[],[],[],phi,[],sigma)

