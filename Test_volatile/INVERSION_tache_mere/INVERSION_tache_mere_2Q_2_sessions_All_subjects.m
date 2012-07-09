% This scripts loads data from the TACHE MERE done by stefano Palminteri
% It then performs the inversion of the DATA using a model in which
% probabilities rather than expected outcome are learned

clear all
close all

%----------- DO THE INVERSION FOR ALL SUBJECTS

controls1=[1001:1015 1103 1106 1108 1109]; % controles agés (Stefano, tache mêre)
INVERSION = cell(1,length(controls1));
data = LOAD_and_FORMAT_data_tache_mere(); % This function loads and formats data for inversion


%controls1=[1002]; % controles agés (Stefano, tache mêre)

controls1 = controls1(1);
i_s =0;
for i_subject = 2
i_s = i_s+1;
i_s = 3;

Y = data{i_s}.data.Y_QL;
U = data{i_s}.data.U_QL;

%%
%----------- SECOND : Performing the inversion with probability learning

Nsessions = 4;
Ntrials = size(Y,2);
dim_output = 1; % the delta reaction-time
dim_data = 2; % reward/index of sequence 
dim = struct('n',2*Nsessions,...  %( 2 (special for RL) * 2 (Qvalues) * Nsessions
             'n_ps', 2,... % number of hidden state per session
             'p',Nsessions*dim_output,... % total output dimension
             'p_ps',dim_output,... % output dimension per sessions
             'n_theta',1,... % evolution parameters
             'n_phi', 1,... % observation parameters
             'n_data_ps',dim_data,... % data dimension per session
             'n_t',Ntrials,...
             'n_sess', Nsessions); %

% Defining Priors
% Priors on parameters (mean and Covariance matrix)
priors.muPhi = zeros(dim.n_phi,1); 
priors.muTheta = zeros(dim.n_theta,1);
priors.SigmaPhi = 1e2*eye(dim.n_phi);
priors.SigmaTheta = 1e2*eye(dim.n_theta);

% Priors on initial 
priors.muX0 = ones(dim.n,1)*0;
priors.SigmaX0 = 0e4*eye(dim.n);
% No state noise for deterministic update rules
priors.a_alpha = Inf;
priors.b_alpha = 0;

% Defining parameters for sessions
in_sessions = struct();
in_sessions.f_fname = @f_Qlearn_2Q;
in_sessions.g_fname = @g_softmax_2Q;
in_sessions.dim = dim;
in_sessions.ind.theta =  ones(Nsessions,1)*[1];
in_sessions.ind.phi =   ones(Nsessions,1)*[1];
in_sessions.binomial = 1;
in_sessions.inF.param_transform.type= 'sigmoid';
in_sessions.inG.param_transform.type= 'exponential';

IsYout = 0 * Y;
[posterior,out] = inversion_multiple_sessions(in_sessions, Y, U, IsYout, priors);


% Saving data after inversion
clear results;
results.posterior = posterior;
results.out = out;
results.in_sessions = in_sessions;
results.Y = Y;
results.U = U;
results.IsYout = IsYout;
results.priors = priors;
INVERSION{i_s}.results = results;

end



save sujet2Q INVERSION
