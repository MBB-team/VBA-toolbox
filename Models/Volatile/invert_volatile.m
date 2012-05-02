%%%%%%%%%%%% simulation
close all
clear all
% Initial states

N =400;



%-------------------------- parameters
% - P: the perceptual model parameters vector, ie. P = [ka;om;th], using
% s2 = ka*x3+om
% s3 = th


ka = 1;
om = 10;
th = 1;

in.thub = 10; % upper bound on theta (lower being 0)
in.kaub = 10; % upper bound on kappa (lower being 0)

% This parameter does not figure in the initial paper.
% When updated, variance on 3rd level state can get negative
% this happens during a reduction of the variance.
% When this happens, an alternative update is performed through a
% rescaling.
in.rf = 1000;

theta = [sigm(ka/in.kaub,struct('INV',1));om;sigm(th/in.thub,struct('INV',1))];

%-------------------------- initial states

x0 = zeros(7,1);
x0(1)=0;%   mu1 = x(1);
x0(2)=0;%   mu2 = x(2);
x0(3)=10;%   sa2 = x(3);
x0(4)=10;%   mu3 = x(4);
x0(5)=10;%   sa3 = x(5);


%-------------------------- inversion
alpha = Inf;
sigma = Inf;
options.binomial = 1;
options.inG = in;
options.inF = in;




%------------------------- Simulating data
Lb = N/2;
ul = [rand(1,Lb/2)<0.1,rand(1,Lb/2)<0.9];
uc = [rand(1,Lb/4)<0.1,rand(1,Lb/4)<0.9]; uc = [uc,uc];
u = [ul,uc];
u=[u, zeros(1,N-size(u,2))];

y = u;
dim_output = 1; % 
dim_data = 1; % reward/index of sequence 
dim = struct('n',7,...  %( 2 (special for RL) * 2 (Qvalues) * Nsessions
             'p',dim_output,... % total output dimension
             'n_theta',3,... % evolution parameters
             'n_phi', 0,... % observation parameters
             'n_t',N);

% Defining Priors
% Priors on parameters (mean and Covariance matrix)
priors.muPhi = zeros(dim.n_phi,1); 
priors.muTheta = zeros(dim.n_theta,1);
priors.SigmaPhi = 1e4*eye(dim.n_phi);
priors.SigmaTheta = 1e4*eye(dim.n_theta);
% Priors on initial 
priors.muX0 = 1*ones(dim.n,1);
priors.SigmaX0 = 0e2*eye(dim.n);
% No state noise for deterministic update rules
priors.a_alpha = Inf;
priors.b_alpha = 0;
% Defining parameters for sessions

options.inF = in;
options.inG = in;
options.priors = priors;
options.DisplayWin = 1;
options.GnFigs = 0;
options.binomial = 1; % Dealing with binary data
options.isYout = zeros(1,N); % Excluding data points

f_fname = @f_Mathys_binary;
g_fname = @g_Mathys_binary;

 [posterior,out] = VBA_NLStateSpaceModel(y,u,f_fname,g_fname,dim,options);

 
 
%x0(1)=0;%   mu1 = x(1);
%x0(2)=0;%   mu2 = x(2);
%x0(3)=10;%   sa2 = x(3);
%x0(4)=10;%   mu3 = x(4);
%x0(5)=10;%   sa3 = x(5);
% posterior.muX(1,:)
% posterior.muX(3,:)
 
% posterior.muPhi
% posterior.muTheta
 
 