clear all 
close all
clc

% Small example to test possible sources of errors in the inversion
% procedure

% Model : unidimensional massive object linked to a linear spring (fixed to zero, equilibrium length =0, constant = k) 
% and subject to a known external force
% random fluctuations affect both position and speed.
% - hidden states (2*1)
%    - speed
%    - position
% - evolution parameters (2*1)
%    - spring constant/mass
%    - dissipation term/mass
% - observation parameters (1*1)
%    - sigmoid parameter
% - input
%    - perturbation strength/mass (1*1)

% Models dynamic.
% Classical physics, dissipation proportional to speed

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% DATA SIMULATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Ntrials = 50;
u = 2*randn(1,Ntrials) + 2; % external force applied to point
theta = [log(4);log(0.1)]; % recall, dissipation
phi = 0; % transform parameter (sigmoid inverse temperature)
x0 = [0;3]; % initial position (speed, position)

inF.dt = 0.1; % delta-time

f_fname = @f_test_errors;
g_fname = @g_test_errors;

alpha = 1e5; % state noise precision
sigma = 1e5; % measurement noise precision

dim_data = 1; % 1D force
dim = struct('n',2,...  % speed, position
             'p',1,... % % 1D position 
             'n_theta',2,... % recall parameter
             'n_phi', 1,... % observation parameter
             'n_t',Ntrials);

options.DisplayWin = 1;
options.GnFigs = 0;
options.binomial = 0; % Dealing with binary data
options.isYout = zeros(1,Ntrials); % Excluding data points
options.inF = inF;
options.dim = dim;
options.backwardLag = 10;

options.priors.iQx =cell(1,Ntrials);
for i  = 1:Ntrials
options.priors.iQx{i} = [0.1,0;0,1];
end

[y,x,x0,eta,e] = simulateNLSS(Ntrials,f_fname,g_fname,theta,phi,u,alpha,sigma,options,x0);

figure
plot(x')
hold on 
plot(y,'r')

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% MODEL INVERSION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Defining Priors
% Priors on parameters (mean and Covariance matrix)
priors.muPhi = zeros(size(phi)); 
priors.muTheta = zeros(size(theta));
priors.SigmaPhi = 1e0*eye(dim.n_phi);
priors.SigmaTheta = 1e0*eye(dim.n_theta);
% Priors on initial 
priors.muX0 = zeros(size(x0));
priors.SigmaX0 = 1e0*eye(dim.n);
% No state noise for deterministic update rules
% priors.a_alpha = 1;
% priors.b_alpha = 1;
priors.a_alpha = 1;
priors.b_alpha = 1;
% Defining parameters for sessions

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% MODEL INVERSION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% Options for inversion
options.DisplayWin = 1;
options.GnFigs = 0;
options.binomial = 0; % Dealing with binary data
options.priors = priors;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% ADDING ERROR TERM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%--- ERRORS IN OBSERVATION FUNCTION FUNCTION
% ERROR TYPE : INDEX out of bound
% #1 : wrong call of parameter
% #2 : wrong call of data
% #3 : wrong call of hidden states
% ERROR TYPE : WRONG calculation
% #4 : erroneous calculation function
% ERROR TYPE : WRONG output
% #5 : output of function is nan/inf
% #6 : output of function has wrong dimensions
errG = 6;

%--- ERRORS IN EVOLUTION  FUNCTION
% ERROR TYPE : INDEX out of bound
% #1 : wrong call of parameter
% #2 : wrong call of data
% #3 : wrong call of hidden states
% ERROR TYPE : WRONG calculation
% #4 : erroneous calculation function
% ERROR TYPE : WRONG output
% #5 : output of function is nan/inf
% #6 : output of function has wrong dimensions
% ERROR TYPE : WRONG input
% #7 : input is nan
errF = 0;


options.inF.error = errF;
options.inG.error = errG;

disp('')
disp(['Error number : ',num2str([errF,errG])])

[posterior,out] = VBA_NLStateSpaceModel(y,u,f_fname,g_fname,dim,options);

displayResults(posterior,out,y,x,x0,theta,phi,Inf,Inf)
