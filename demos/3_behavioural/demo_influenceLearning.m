% This script simulates and inverts an "influence learner"
% This model is based upon [Hampton et al., 2008].
% Let p=P(o=1) be the agent's prediction of the other's next move, ie the
% probability that the other will pick the first alternative option. The
% "influence learning" rule can be written as follows:
% p <- p + eta*PE1 + lambda*p*(1-p)*PE2
% where o is the other's last move, PE1=o-p is the first-order prediction
% error, PE2 is the second-order prediction error, eta is the weight of
% PE1 and lambda is the weight of PE2.
% Note: PE2 drives the "influence" learning. If set to 0, then the
% influence learner reduces to a simple RL agent (with counterfactuals).

% clear all
close all
clc

% simulation parameters
payoffTable = cat(3,[1,0;0,1],[0,1;1,0]); % game payoff matrix (here: hide-and-seek)
role = 1;  % player's role (here: 1=seeker, 2=hider) 
options.inF = struct('game',payoffTable,'player',role);
options.inG = struct('game',payoffTable,'player',role);
dim.n = 1;
dim.n_theta = 3;
dim.n_phi = 2;



%% simulate sequence of k-ToM choices
x0 = [0]; % log-odds of P(o=1)
theta = [1;1;0]; % weight (PE1), weight (PE2), opponent's temp
phi = [-1;0]; % (log-) temperature, bias
N = 50; % number of trials
o = VBA_random ('Bernoulli', 0.75, 1, N); % opponent's choices (here dummy binomial sampling)
a = NaN(1,N);
gx = NaN(1,N);
x = zeros(dim.n,N+1);
x(:,1) = x0; %initialize hidden states
for i=1:N
    gx(i)= g_Hampton(x(:,i)',phi, [],options.inG) ; 
    a(i)= gx(i)>.5;
    x(:,i+1)= f_Hampton(x(:,i),theta, [o(i);a(i)],options.inF);
end
% figure,plot(x')

%% invert "influence learning" model given sequence of agent's choices
options.skipf = zeros(1,N);
options.skipf(1) = 1;
options.sources = struct('type',1 ,'out', 1); % one binomial observation;
options.priors.SigmaTheta = 1e2*eye(dim.n_theta);
f_fname = @f_Hampton;
g_fname = @g_Hampton;
u = [zeros(2,1),[o;a]];
[posterior,out] = VBA_NLStateSpaceModel(a,u,f_fname,g_fname,dim,options);



