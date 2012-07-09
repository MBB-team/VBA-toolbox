% Inversion of the synthetic data from and with a hyperbolic discount
% model with softmax decision in a binary decision task

close all
clear all
clc

%--------------------- Task
% Two alternatives choice task
% - At each Trial, two alternatives are proposed by the experimenter
%   Each of which has two dimensions (value, time of delivery)
%   Delay between alternatives is kept constant 
% - The Subject chooses one among the two alternatives
%--------------------- Model
% Hyperbolic discount v(d) = 1 / (1 + k*d) + softmax decision rule
% Hidden states : None
% Observed variable : chosen alternative (1*Ntrials)
% Parameters : k, be
%---------------------------



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% INVERSION OF MODEL ONE


g_fname = @g_1Dhyp; % observation function, Parameters : [K, beta] => n_phi = 2
N = 100; % number of trials
dim = struct('n',0,'n_theta',0,'n_phi',2,'p',N,'n_t',1);
% Priors on parameters (mean and Covariance matrix)
priors.muPhi = zeros(dim.n_phi,1); 
priors.SigmaPhi = 1e4*eye(dim.n_phi);
%priors.SigmaPhi(end,end) = 0; % Do not infer beta!
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

%----------------- Simulating data
% simulation parameters
K = 1/15; % discount parameter
beta = 0.5; % inverse temperature
% Experimenter data
T = zeros(2,N); % time of reception of alternatives
T(1,:) = zeros(1,N); % random choice (uniform)  
T(2,:) = T(1,:) + rand(1,N)*30; % T2 = T1 + random delay  
OV = zeros(2,N); % value of alternatives
OV(1,:) = 10; % fixed objective value
OV(2,:) = 20; % fixed objective value
% simulation
inG.T = T; % 2 * N (times of reception)
inG.V = OV;  % 2 * N  (objective values)
options.inG = inG; 
[y,x,x0,eta,e] = simulateNLSS(N,[],g_fname,[],[log(K);log(beta)],[],Inf,Inf,options,[]);

%-----------------Model inversion
[posterior,out] = VBA_NLStateSpaceModel(y,[],[],g_fname,dim,options);  % Inversion function
 phi = [log(K),log(beta)]; 
displayResults(posterior,out,y,[],[],[],phi,Inf,Inf)


disp('//////////  Obtained parameters for hyperbolic model')
disp([' K = ', num2str(exp(posterior.muPhi(1)))])
disp([' beta = ', num2str(exp(posterior.muPhi(2)))])
    
    
