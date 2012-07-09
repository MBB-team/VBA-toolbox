clear all
close all
clc 

P = [sigm(0.2,struct('INV',1)),0,1,10,1];

N = 500;
u = [];
% generate from mixture
for i = 1 : N
    u(1,i) = g_MoG_1D([],P,[],[]);
end

%%

f_fname = @g_classify_MoG_1D_from_classification; 

g_fname = @g_classify_MoG_1D_from_classification; 
dim = struct('n',0,...
             'n_theta',0,...
             'n_phi',5,...
             'p',1,...
             'n_t',N);


% Priors on parameters (mean and Covariance matrix)
priors.muPhi = zeros(dim.n_phi,1); 
priors.SigmaPhi = 1e4*eye(dim.n_phi);

priors.SigmaPhi(1,1) = 1;
priors.SigmaPhi(3,3) = 0;
priors.SigmaPhi(5,5) = 0;

% No state noise for deterministic update rules
priors.a_alpha = Inf;
priors.b_alpha = 0;
% Options for inversion
options.priors = priors;
options.DisplayWin = 1;
options.GnFigs = 0;
options.binomial = 1; % Dealing with binary data
options.dim = dim;
options.verbose = 0;
% Experimenter parameters


% Simulating
phi = [sigm(0.5,struct('INV',1)),0,1,10,1];
[y,x,x0,eta,e] = simulateNLSS(dim.n_t,[],g_fname,[],phi,u,Inf,Inf,options,[]);



[posterior,out] = VBA_NLStateSpaceModel(y,u,[],g_fname,dim,options);  % Inversion function
 
%%
%%
clear all
close all

P = [sigm(0.2,struct('INV',1)),0,1,10,1];

N = 100;
u = [];
% generate from mixture
for i = 1 : N
    u(1,i) = g_MoG_1D([],P,[],[]);
end
figure
hist(u,20)

%%
g_fname = @g_MoG_1D; 
dim = struct('n',0,...
             'n_theta',0,...
             'n_phi',5,...
             'p',1,...
             'n_t',N);


% Priors on parameters (mean and Covariance matrix)
priors.muPhi = zeros(dim.n_phi,1); 
priors.muPhi =[sigm(0.5,struct('INV',1));0;1;10;1];

priors.SigmaPhi = 1e4*eye(dim.n_phi);



% No state noise for deterministic update rules
priors.a_alpha = Inf;
priors.b_alpha = 0;
% Options for inversion
options.priors = priors;
options.DisplayWin = 1;
options.GnFigs = 0;
options.binomial = 0; % Dealing with CONTINUOUS data
options.dim = dim;
options.verbose = 0;
% Experimenter parameters


% Simulating
phi = [sigm(0.5,struct('INV',1)),0,1,10,1];
[y,x,x0,eta,e] = simulateNLSS(dim.n_t,[],g_fname,[],phi,u,Inf,Inf,options,[]);

hist(y,20)

%%

[posterior,out] = VBA_NLStateSpaceModel(y,u,[],g_fname,dim,options);  % Inversion function

displayResults(posterior,out,y,x,x0,[],phi,Inf,Inf)

