% This script is a test of the function to assess multiple parametric
% modulations

% In this example, I consider
% - a RL model of two parameters
% - 2 inputs that I will test as modulators of the two parameters

close all 
clear all

% ---- Data simulation ----
Ntrials = 200;

% -- Simulate input
u = [zeros(2,Ntrials);... % first two lines left blank for now  
    randn(1,Ntrials); % input 1
    randn(1,Ntrials); % input 2
    ];
% u : 
% - 1 : action
% - 2 : reward
% - 3,4 : modulators

% -- Define Contingencies

R = rand(2,Ntrials); % Rewards for each alternatives
R(1,1:Ntrials/2) = R(1,1:Ntrials/2)<0.7;R(1,Ntrials/2+1:end) = R(1,Ntrials/2+1:end)<0.3;
R(2,1:Ntrials/2) = R(2,1:Ntrials/2)<0.3;R(2,Ntrials/2+1:end) = R(2,Ntrials/2+1:end)<0.7;

% -- Simulate data
% RL model: 
% - learning rate alpha = 0.3, modulation by u(1) with linear parameter
% 0.1;
% - decision noise beta = 1, modulation by u(2) with linear paramter 0.1

x0 = [0;0];
n_t = Ntrials;
fb.inH.u0 = R;
fb.h_fname = @h_reward_2Q;
fb.indy = 1;
fb.indfb = [2]; % indices where put feedbacks in the experimenter data matrix u
f_fname = @f_Qlearn_2Q;
g_fname = @g_softmax_2Q;

alpha =0.2; theta = [sigm(alpha,struct('INV',1))];
beta = 1; phi = [log(beta)];

alpha = Inf;
sigma =  Inf;
dim_output = 1; % the delta reaction-time
dim_data = 3; % index of sequence(not used for sim)/ rewards for both alternatives
dim = struct('n',2,...  %( 2 (Qvalues) * Nsessions
             'p',dim_output,... % total output dimension
             'n_theta',1,... % evolution parameters
             'n_phi', 1,... % observation parameters
             'n_t',Ntrials);

options.inF.param_transform.type = 'modified sigmoid';
options.inF.param_transform.a = 0;
options.inF.param_transform.b = 1;
options.inG.param_transform.type = 'modified sigmoid';
options.inG.param_transform.a = 1;
options.inG.param_transform.b = 5;

options.DisplayWin = 1;
options.GnFigs = 0;
options.binomial = 1; % Dealing with binary data
options.isYout = zeros(1,Ntrials); % Excluding data points
options.dim = dim;
              
options.skipf = zeros(1,n_t);
options.skipf(1) = 1; % apply identity mapping from x0 to x1.

[y,x,x0,eta,e,u] = simulateNLSS_fb(n_t,f_fname,g_fname,theta,phi,u,Inf,Inf,options,x0,fb);

plot(x')
hold on
plot(y,'x')

%%



phi = cell(1,dim.n_theta);

% Declaring modulation parameters

mod.phi = cell(1,1);
mod.phi{1}.indp = 1; % indice of parameter
mod.phi{1}.indu = [2,3]; % indice of modulatory inputs

mod.theta = cell(1,1);
mod.theta{1}.indp = 1; % indice of parameter
mod.theta{1}.indu = [2,3]; % indice of modulatory inputs

% Building extended model for modulation

[ f_fname_m,g_fname_m,dim_m,options_m, mod ] = build_modulated_model( f_fname,g_fname,dim,options,mod );
%%

% -- declaring default priors

% Priors on parameters (mean and Covariance matrix)
priors.muPhi = zeros(dim_m.n_phi,1); 
priors.muTheta = zeros(dim_m.n_theta,1);
priors.SigmaPhi = 1e0*eye(dim_m.n_phi);
priors.SigmaTheta = 1e0*eye(dim_m.n_theta);
% Priors on initial 
priors.muX0 = ones(dim_m.n,1)*0;
priors.SigmaX0 = 0e4*eye(dim_m.n);
% No state noise for deterministic update rules
priors.a_alpha = Inf;
priors.b_alpha = 0;

options_m.DisplayWin = 1;
options_m.GnFigs = 0;
options_m.binomial = 1; % Dealing with binary data
options_m.priors = priors;



%%
[inversions,inv_order,p_m,LogEv] = VBA_Inversion_modulation2(y,u,f_fname_m,g_fname_m,dim_m,options_m,mod);

