% This scripts loads data from the TACHE MERE done by stefano Palminteri
% It then performs the inversion of the DATA using a model in which
% probabilities rather than expected outcome are learned

clear all
close all

%----------- FIRST : Loading the data & formating for inversion

root = 'C:\Users\vincent.adam\Desktop\Vincent ADAM\Matlab\Demo_toolbox\Data\Data_tache_mere\VINCENT\'; % location of data folder
subfolder = 'CONTROLS_1\';

i_subject = 1002; % index of subject


%-- LOADING THE TWO SESSIONS --%

SESSIONS = [1,3];
Y = [];
U = [];

for i_session = SESSIONS; % index of session

filename = ['ILtestSub',num2str(i_subject),'Session',num2str(i_session)]; % as saved by stefano
file_root = [root, subfolder, filename]; % position of file on hard-drive

load(file_root); % loading file

%-----NOTES on how file is organized
% session1 = training session
% session2-3 = test sessions
%-- DATA
% column 1 = index of session
% column 2 = index of trial
% column 3 = condition (1=Gain (1€), 2=Neutral (0€), 3=Loss (-1€))
% column 4 = checktime
% column 5 = 1=Go, -1=NoGo 
% column 6 = choice (1=correct, -1=incorrect)
% (correct = most rewarded option in cond=1 | less punished option in cond=3)
% colonne 7 = feedback (cond 1 : 1 = 1€ , -1=0€ | cond 3, 1=0€, -1=-1€)
% colonne 8 = reaction time, (millisec if Go, 0 if NoGo)
% colonne 9 = number of hesitations

%---- Formating data for further inversion

%-----NOTES on how data should be formated
% y: binary choice 
% u: 
%   - subject choice (binary)
%   - outcome category (binary)
% 
% outcomes values must be specified in options:
%   - inG.o1;
%   - inG.o2;


%=== second, let's consider the two conditions together : gains and loss.
cond = data(:,3); % vector of conditions
i_cond1 = find(cond==1); % indices of trial where condition is 1

% outputs
y = data(i_cond1,6)'>0; % binary choice (0: correct, 1:incorrect)
% experimenter data
u = [y;... % binary choice (0: correct, 1:incorrect)
     data(i_cond1,7)'<0;... % binary outcome (0: best, 1:worst)
     data(i_cond1,7)'>0]; 
 
 % REMARK : Contrary to the probability learning model with known outcomes, reinforcement learning needs to be fed by 
 % - choice (binary)
 % - outcomes of both choices (this was chosen to enable simulations of the
 % model

Y = [Y;y]; 
U = [U;u]; 

end
%%
%----------- SECOND : Performing the inversion with probability learning

Nsessions = 2;
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
in_sessions.inF = [];
in_sessions.inG = [];

IsYout = 0 * Y;
[posterior,out] = inversion_multiple_sessions(in_sessions, Y, U, IsYout, priors);
