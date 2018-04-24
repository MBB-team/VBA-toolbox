% demo for stochastic system inversion given binomial data

clear variables
close all


p = 5e2;
u = 1e0*randn(p,1);
% u = sort(u);

% get parameters trajectory and binomial sampling probabilities over time
g_fname = @g_sigm_binomial;
f_fname = @f_ARn;
theta = [-5;-6;3;4];
x0 = [1;-4];
options.inG.x = 1;
[sx,x,x0,eta,e] = simulateNLSS(p,f_fname,g_fname,theta,[],u(:)',Inf,Inf,options,x0);

% sample binomial data
y = zeros(p,1);
for t=1:p
    [y(t)] = sampleFromArbitraryP([sx(t),1-sx(t)]',[1,0]',1);
end

figure
subplot(2,1,1)
plot(x')
title('sigmoid parameters over time (trials)')
legend({'sigmoid slope','sigmoid inflextion point'})
subplot(2,1,2)
plot(sx,'r')
hold on
plot(y,'k.')
title('choices and their likelihood over time (trials)')
legend({'sampling proba: p(y|x)','data samples: y'})


dim.n_phi = 0;                  % nb of observation parameters
dim.n_theta = 4;                % nb of evolution parameters
dim.n=2;                        % nb of hidden states

% Call inversion routine
options.sources = struct('type',1 ,'out', 1); % one binomial observation;
priors.muTheta = 0.*ones(dim.n_theta,1);
priors.SigmaTheta = 1e2*eye(dim.n_theta);
priors.SigmaX0 = 1e2*eye(dim.n);
% priors.a_alpha = 1;
% priors.b_alpha = 1;
options.priors = priors;
options.DisplayWin = 1;
options.GnFigs = 0;
options.backwardLag = 16;
options.checkGrads = 0;

[posterior,out] = VBA_NLStateSpaceModel(y(:)',0.*u(:)',f_fname,g_fname,dim,options);

displayResults(posterior,out,y,x,x0,theta,[],[],[]);
