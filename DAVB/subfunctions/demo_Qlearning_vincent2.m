
%%%%%%%%%%%%% Description du modèle.

% Le sujet doit choisir parmi 2 actions,
% A chaque action est attribuée une valeur, que le sujet apprend au travers
% d'un modèle de rescorla wagner (learning rate alpha fixe)
% A chaque essai, le sujet choisit une des deux actions. Le choix est
% effectué par une règle de décision softmax (gouvernée par un seul
% paramètre beta = température inverse)


% hidden states : Q-values (2*Ntrials)
% observed variable : chosen action (1*Ntrials)
% parameters : alpha,beta. (learning rate and inverse temperature)

close all
clear variables
clc


% simulation parameters
alpha = 0.3;
beta = 4;
u = [ones(1,50);zeros(1,50)]; % possible feedbacks
u = [u,~u,u,~u,u,~u]; % with inversions
N = length(u);

% pre-allocate variables
Q = zeros(2,N+1); % Q-values
Q(:,1) = 0.5;
A = zeros(1,N); % chosen action time series
p = zeros(2,N); % action likelihood
r = zeros(1,N); % obtained feedback

for t = 1 : N
    
    % action emission
    p1 = 1/(1+exp(beta*(Q(2,t)-Q(1,t))));
    p(:,t) = [p1;1-p1];
    a = sampleFromArbitraryP(p(:,t)',[1,0],1);   
    A(t)=a;
    % feedback
    r(t) = u(a+1,t);
    % Q-values update
    Q(:,t+1)= Q(:,t);
    Q(a+1,t+1)=Q(a+1,t)+alpha*(r(t)-Q(a+1,t));
    
end

% display simulated data
figure

subplot(2,1,1)
hold on
plot(1:N,u(1,:),'r')
plot(1:N,u(2,:),'b')
plot(1:N,Q(1,1:end-1),'.-r')
plot(1:N,Q(2,1:end-1),'.-b')
title('Q-values')

subplot(2,1,2)
hold on
plot(p(1,:),'r')
plot(p(2,:),'b')
plot(1:N,A,'xk')
title('p(action|P) and emitted action')




u = [A;r]; % Input (2*Ntrials)
y = A; % Output : chosen action

f_fname = @f_Qlearn; % P_f = alpha => n_theta = 1
g_fname = @g_softmax; % P_g = beta => n_phi = 1

dim = struct('n',2,'n_theta',1,'n_phi',1);

priors.muPhi = zeros(dim.n_phi,1);
priors.muTheta = zeros(dim.n_theta,1);
priors.muX0 = zeros(2,1);
priors.SigmaPhi = 1e4*eye(dim.n_phi);
priors.SigmaTheta = 1e4*eye(dim.n_theta);
priors.SigmaX0 = 1e4*eye(dim.n);
priors.a_alpha = Inf;
priors.b_alpha = 0;

options.priors = priors;
options.DisplayWin = 1;
options.GnFigs = 0;
options.binomial = 1;

% use simulateNLSS.m, but note:
% - parameter alpha is defined as sigm(theta(1)) (c.f. @f_Qlearn)
% - parameter beta is defined as exp(phi(1)) (c.f. @g_softmax)
theta = sigm(alpha,struct('INV',1)); 
phi = log(beta); 
Q0 = [0.5;0.5]; % Initialisation des Qvalues
[pc,Q,x0,eta,e] = simulateNLSS(...
    N,f_fname,g_fname,theta,phi,u,Inf,Inf,options,Q0);
%%%
figure
plot(pc-e,'r')
hold on
plot(y,'kx')
legend({'p(y=1|theta,phi,m)','binomial data samples'})
getSubplots
pause


options.isYout = zeros(1,size(y,2));
options.isYout(2) = 1;

[posterior,out] = VBA_NLStateSpaceModel(y,u,f_fname,g_fname,dim,options)

displayResults(posterior,out,y,Q,Q0,theta,phi,Inf,Inf)


