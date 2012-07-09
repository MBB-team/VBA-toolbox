% This scripts loads data from the TACHE MERE done by stefano Palminteri
% It then performs the inversion of the DATA using a model in which
% probabilities rather than expected outcome are learned

clear all
close all

%----------- FIRST : Loading the data & formating for inversion

root = 'C:\Users\vincent.adam\Desktop\Vincent ADAM\Matlab\Demo_toolbox\Data\Data_tache_mere\VINCENT\'; % location of data folder
subfolder = 'CONTROLS_1\';

i_subject = 1002; % index of subject
i_session = 3; % index of session

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


%=== first, let's consider the gain condition alone, that is cond 1.
cond = data(:,3); % vector of conditions
i_cond1 = find(cond==1); % indices of trial where condition is 1

% outputs
y = data(i_cond1,6)'>0; % binary choice (0: correct, 1:incorrect)
% experimenter data
u = [y;... % binary choice (0: correct, 1:incorrect)
     data(i_cond1,7)'>0]; % binary outcome (0: best, 1:worst)

%----------- SECOND : Performing the inversion with probability learning


% evolution, observation and feedback functions
f_fname = @f_OpLearn_2p; % evolution function for the learning of probability
g_fname = @g_softmax_EU_2p; % softmax decision based on probabilities of outcome 1
% parameters : 1=inverse temperature, 2= asymmetry parameter in prospect
% theory
h_fname = @h_choice_outcome_2p; % feedback is the outcome of the chosen action

% allocate feedback struture for simulations

% defining the utility function
u_fname = @u_prospect_theory; % handle of utility function
inG.u_fname = u_fname;
inG.o1 = 1;
inG.o2 = 0;

% simulation parameters % See  Mathys, Daunizeau et al. 2010 for a
% detailled description of parameters
inF.lev2 = 1; % remove 3rd level (volatility learning)
inF.kaub = 1.4;
inF.thub = 1;
inF.rf = -1;
inG.respmod = 'taylor';

% choose initial conditions
x0 = repmat([0.5;0;0;1;log(4)],2,1);

dim = struct('n',2*5,... % number of hidden states in the probability learning model
             'n_theta',3,...
             'n_phi',2);

priors.muPhi = zeros(dim.n_phi,1);
priors.muTheta = [0;-4;0];
priors.muX0 = x0;
priors.SigmaPhi = 1e2*eye(dim.n_phi);
priors.SigmaTheta = 1e2*eye(dim.n_theta);
priors.SigmaX0 = 0e1*eye(dim.n);
priors.a_alpha = Inf;
priors.b_alpha = 0;

options.priors = priors;
options.binomial = 1;
options.inF = inF;
options.inG = inG;
options.skipf = zeros(1,length(u));
options.skipf(1) = 1; % apply identity mapping from x0 to x1.
options.isYout = zeros(1,length(y));

[posterior,out] = VBA_NLStateSpaceModel(y,u,f_fname,g_fname,dim,options);
