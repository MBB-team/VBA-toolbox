clear all;
close all;
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% Inverting data of all groups together
%%%%%%%%%%%%%%% Separate loss and gains
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[Y,U,IsYout] = Simulate_data_Pessiglione2006();
Nsessions_tot = size(Y,1);
 Ntrials =  size(Y,2);

dim_data = 4;
dim_output = 1;
         
dim = struct('n',4*Nsessions_tot,...
             'n_ps', 4,... % number of hidden state per session
             'p',Nsessions_tot*dim_output,... % total output dimension
             'p_ps',dim_output,... % output dimension per sessions
             'n_theta',2,... % evolution parameters
             'n_phi', 2,... % observation parameters
             'n_data_ps',dim_data,... % data dimension per session
             'n_t',Ntrials,...
             'n_sess', Nsessions_tot); %

options.binomial = 1; % Dealing with binary data

% Defining Priors
% Priors on parameters (mean and Covariance matrix)
priors.muPhi = zeros(dim.n_phi,1); 
priors.muTheta = zeros(dim.n_theta,1);
priors.SigmaPhi = 1e4*eye(dim.n_phi);
priors.SigmaTheta = 1e4*eye(dim.n_theta);
% Priors on initial 
priors.muX0 = ones(dim.n,1)*0;
priors.SigmaX0 = 0e1*eye(dim.n);

% No state noise for deterministic update rules

priors.a_alpha = Inf;
priors.b_alpha = 0;

% in_session
in_sessions = struct();
in_sessions.f_fname = @f_Qlearn_Pessiglione2006_2Q_M3;
in_sessions.g_fname = @g_softmax_Pessiglione2006_2Q_M3;
in_sessions.dim = dim;
in_sessions.binomial = 1;

in_sessions.ind.theta =  ones(Nsessions_tot,1)*[1:2]; % same theta for each session
in_sessions.ind.phi =  ones(Nsessions_tot,1)*[1:2];   % same phi for each session

% Performing the inversion
[posterior,out] = inversion_multiple_sessions(in_sessions, ~Y, U, IsYout, priors);
%  pause 





