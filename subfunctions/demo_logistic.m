% this function exmeplifies the use of the logistic model

clear variables
close all

dim.p = 128;
dim.n = 0;
dim.n_t = 1;
dim.n_theta = 0;
dim.n_phi = 4;
inG.X = [ones(dim.p,1),randn(dim.p,dim.n_phi-1)];

% Choose basic settings for simulations
phi = ones(dim.n_phi,1);
g_fname = @g_logistic;
options.inG = inG;
options.sources = struct('type',1 ,'out', 1); % one binomial observation;
options.dim = dim;
options.checkGrads = 0;
% Build simulated observations 
[y,x,x0,eta,e] = simulateNLSS(dim.n_t,[],g_fname,[],phi,[],[],[],options);
% display time series of hidden states and observations
figure,
plot(y','ro')
hold on;
plot((y-e)')


% Call inversion routine
[posterior,out] = VBA_NLStateSpaceModel(y,[],[],g_fname,dim,options);
% Display results
displayResults(posterior,out,y,[],[],[],phi,[],[]);