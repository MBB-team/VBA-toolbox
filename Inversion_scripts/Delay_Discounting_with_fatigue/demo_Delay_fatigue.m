% binomial data inversion of two models for delay discounting
%%% This script is divided into the following parts : 
% - Simulation of data from one of the two models
% - Inversion of the generated data for both model
% - Model comparison based on the 

%SCRIPT POUR L'EXEMPLE DE LA DOC/


close all
clear all
clc

Nmodel = 3;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DESCRIPTION OF THE TASK

% Two alternatives choice task
% - At each Trial, two alternatives are proposed by the experimenter
%   Each of which has two dimensions (value, time of delivery)
%   Delay between alternatives is kept constant 
% - The Subject chooses one among the two alternatives
% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DESCRIPTION OF MODEL ONE
% Hyperbolic discount v(d) = 1 / (1 + k*d) + softmax decision rule
% Hidden states : None
% Observed variable : chosen alternative (1*Ntrials)
% Parameters : k, be

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SIMULATION OF DATA WITH MODEL ONE


% simulation parameters
beta = 3;


% Experimenter data
N = 100;

T = zeros(2,N); % time of reception of alternatives
T(1,:) = zeros(1,N); % random choice (uniform)  
T(2,:) = 40+floor((0.5-rand(1,N))*10); % T2 = T1 + delay  
OV = zeros(2,N); % value of alternatives
OV(1,:) = 10; % random choice (uniform)
OV(2,:) = 20; % random choice of values

% Simulating data
K1 = 0.02;
Kend = 0.03;
K = (Kend-K1)/(N-1).*[1:N] + (K1*N-Kend)/(N-1);


SV = (OV)./(1+[K;K].*T);
p1 = 1./(1+exp(-beta*(SV(1,:)-SV(2,:))));
a = rand(1,N)<p1; % subject choices



%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% INVERSION OF MODEL ONE : No fatigue


g_fname = @g_1Dhyp; % observation function, Parameters : [K, beta] => n_phi = 2
dim = struct('n',0,'n_theta',0,'n_phi',2,'p',N,'n_t',1);
% n : hidden states
% n_theta : number of parameters of the evolution function (no evolution
% here)
% n_phi : number of parameters of the observation function

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
% Experimenter parameters
inG.T = T; % 2 * N (times of reception)
inG.V = OV;  % 2 * N  (objective values)
options.inG = inG; 



y = a'; % subject choices

[posterior1,out1] = VBA_NLStateSpaceModel(y,[],[],g_fname,dim,options);  % Inversion function
 phi = [log(K),log(beta)]; 
displayResults(posterior1,out1,y,[],[],[],phi,Inf,Inf)


disp('//////////  Obtained parameters for hyperbolic model')
disp([' K = ', num2str(posterior1.muPhi(1))])
disp([' beta = ', num2str(exp(posterior1.muPhi(2)))])
    
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% INVERSION OF MODEL TWO : Fatigue on discount parameter K


g_fname = @g_1Dhyp_fatigueK; % observation function, Parameters : [K, beta] => n_phi = 2
dim = struct('n',0,'n_theta',0,'n_phi',3,'p',N,'n_t',1);
% n : hidden states
% n_theta : number of parameters of the evolution function (no evolution
% here)
% n_phi : number of parameters of the observation function

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
% Experimenter parameters
inG.T = T; % 2 * N (times of reception)
inG.V = OV;  % 2 * N  (objective values)
inG.N = N;

options.inG = inG; 



y = a'; % subject choices

[posterior2,out2] = VBA_NLStateSpaceModel(y,[],[],g_fname,dim,options);  % Inversion function
 phi = [log(K1),log(Kend),log(beta)]; 
displayResults(posterior2,out2,y,[],[],[],phi,Inf,Inf)


disp('//////////  Obtained parameters for hyperbolic model')
disp([' K = ', num2str(exp(posterior2.muPhi(1:2)'))])
disp([' beta = ', num2str(exp(posterior2.muPhi(3)))])
    
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% INVERSION OF MODEL TWO : Fatigue on inverse temperature B


g_fname = @g_1Dhyp_fatigueB; % observation function, Parameters : [K, beta] => n_phi = 2
dim = struct('n',0,'n_theta',0,'n_phi',3,'p',N,'n_t',1);
% n : hidden states
% n_theta : number of parameters of the evolution function (no evolution
% here)
% n_phi : number of parameters of the observation function

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
% Experimenter parameters
inG.T = T; % 2 * N (times of reception)
inG.V = OV;  % 2 * N  (objective values)
inG.N = N;

options.inG = inG; 



y = a'; % subject choices

[posterior3,out3] = VBA_NLStateSpaceModel(y,[],[],g_fname,dim,options);  % Inversion function
 phi = [log(K1),log(beta),log(beta)]; 
displayResults(posterior3,out3,y,[],[],[],phi,Inf,Inf)

    
    figure
   bar( [out1.F,out2.F,out3.F])
xlabel('models')
title('Log evidence')


