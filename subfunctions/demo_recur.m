% this script demonstrates the simulation and inversion of a k-ToM learner

clear all
close all
clc

K = 2; % depth of k-ToM's recursive beliefs
payoffTable = cat(3,[1,0;0,1],[0,1;1,0]); % game payoff matrix (here: hide-and-seek)
role = 1;  % player's role (here: 1=seeker, 2=hider)
diluteP = 1; % partial forgetting of opponent's level (only for 2-ToM and above)?
[options,dim] = prepare_kToM(K,payoffTable,role,diluteP);


%% simulate sequence of k-ToM choices
phi = [-2;0]; % temperature and bias
theta = -1; % prior volatility
if diluteP==1
    theta = [theta;-5];
end
N = 250; % number of trials
o = bernoulli(.5,N)'; % opponent's choices (here dummy binomial sampling)
tic
a = NaN(1,N); % agent's choices
gx = NaN(1,N);
x = zeros(dim.n,N+1);
[x(:,1)] = f_kToM(options.priors.muX0,theta,[],options.inF); %initialize hidden states
for i=1:N
    gx(i)= g_kToM(x(:,i),phi, [],options.inG) ; 
    a(i)= gx(i)>.5;
    x(:,i+1)= f_kToM(x(:,i),theta, [o(i);a(i)],options.inF);
end
toc
figure,subplot(2,1,1),plot(x')
subplot(2,1,2),plot(a,'b.')
hold on
plot([1,N],[0.5,0.5],'k--')
plot(gx,'r')
unwrapKTOM(x,options.inG)
return

options.skipf = zeros(1,N);
options.skipf(1) = 1;
options.binomial = 1;
options.DisplayWin = 1;
options.priors.SigmaTheta = 1e2*eye(dim.n_theta); % relax evol param
f_fname = @f_kToM;
g_fname = @g_kToM;
u = [zeros(2,1),[o;a]];
[posterior,out] = VBA_NLStateSpaceModel(a,u,f_fname,g_fname,dim,options);



