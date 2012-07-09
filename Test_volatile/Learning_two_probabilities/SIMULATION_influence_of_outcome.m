%----------------------------------------------------------------------
%
% SIMULATING SCRIPT
% INSTRUMENTAL LEARNING TASK
% 1 - Learning probabilities
% 2 - Using a utility function including loss aversion
% 3 - Sotmax decision on expected utility
%
%----------------------------------------------------------------------
close all
clear variables
clc

%-----------------------------------------------------------------------
% PARAMETERS OF SIMULATION

% Probabilities
% Softmax decision inverse temperature parameter

 p1_1 = 0.7; % probability of outcome 1 for choice 1
 p1_2 = 0.3; % probability of outcome 1 for choice 2
 o1 = 1; % value of outcome 1
 o2 = 0; % value of outcome 2  
 phi = [log(2),... inverse temperature
       1];       % Asymetry in utility of prospect theory
 theta = [1;-4;-1];

%-----------------------------------------------------------------------
Ntrials = 100;


O1 = [-5,-1,0,1,5];
LC = zeros(length(O1),Ntrials);

%------- For different values of outcome for the first option (the second
%being kept equal to zero

i_o1 = 0;
for o1 = O1
i_o1 = i_o1 + 1;

Nit = 50;
Y = zeros(Nit,Ntrials);


%------ ITERATING simulation


for it = 1:Nit


% In this task : two actions lead to the same outcomes (binary) with
% different probabilities.  

% 1. Model samples according to a decision rule (see 2.) and track the outcome probabilities of both alternatives 
% 2. Model selects an action based on their probability of beign rewarded :
% Softmax decision on outcome probabilities


% evolution, observation and feedback functions
f_fname = @f_OpLearn_2p; % evolution function for the learning of probability
g_fname = @g_softmax_EU_2p; % softmax decision based on probabilities of outcome 1
% parameters : 1=inverse temperature, 2= asymmetry parameter in prospect
% theory
h_fname = @h_choice_outcome_2p; % feedback is the outcome of the chosen action

% allocate feedback struture for simulations
u1 = [rand(1,Ntrials)<p1_1];  % outcome for alternative 1;
u2 = [rand(1,Ntrials)<p1_2];  % outcome for alternative 2;
fb.inH.u0 = [u1;u2]; % definition of the binary time-series to be predicted
fb.h_fname = h_fname;
fb.indy = 1; % where to write subject choice in vector u
fb.indfb = 2; % where to write subject choice 

% defining the utility function
u_fname = @u_prospect_theory; % handle of utility function
inG.u_fname = u_fname;
inG.o1 = o1;
inG.o2 = o2;

% simulation parameters % See  Mathys, Daunizeau et al. 2010 for a
% detailled description of parameters
inF.lev2 = 1; % remove 3rd level (volatility learning)
inF.kaub = 1.4;
inF.thub = 1;
inF.rf = -1;
inG.respmod = 'taylor';

% choose initial conditions
x0 = repmat([0.5;0;0;1;log(4)],2,1);
u = zeros(2,size(fb.inH.u0,2)+1);

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


% Simulate model and format data for further inversion
[y,x,x0,eta,e,u] = simulateNLSS_fb(length(u),f_fname,g_fname,theta,phi,u,Inf,Inf,options,x0,fb);

Y(it,:) = y(1:end-1);

end



LC(i_o1,:) = mean(Y);

end

figure
hold on
col = 'rgbcymk';

for i = 1 : length(O1)   
    plot(LC(i,:),col(i))
end 
legend('-5','-1','0','1','5')


%%
figure
hold on
plot(y-e,'r')
plot(y,'kx')

plot(sgm(x(2,:),1),'--b')
plot(sgm(x(7,:),1),'--k')
legend({'p(y=1|theta,phi,m)','binomial data samples','p(o=1|a=1)','p(o=1|a=2)'})

getSubplots


figure
hold on
plot(x')
title('hidden states')
% 
% 
% pause
% 
% % Model inversion
% options.isYout = zeros(1,size(y,2));
% [posterior,out] = VBA_NLStateSpaceModel(y,u,f_fname,g_fname,dim,options);
% displayResults(posterior,out,y,x,x0,theta,phi,Inf,Inf)
% 
% 
