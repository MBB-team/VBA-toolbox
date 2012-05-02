% demo for binomial data inversion with adaptation (learning) effect

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
addpath('../sampling')
y = zeros(p,1);
seed = 1e4*rand;
for t=1:p
    try
        [y(t),seed] = binomial_sample(1,sx(t),seed);
    catch
        [y(t)] = sampleFromArbitraryP([sx(t),1-sx(t)],[1,0],1);
    end
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
options.binomial = 1;
% options.gradF = 1;
priors.muTheta = zeros(dim.n_theta,1);
priors.SigmaTheta = 1e4*eye(dim.n_theta);
priors.SigmaX0 = 1e4*eye(dim.n);
priors.a_alpha = Inf;
priors.b_alpha = 0;
options.priors = priors;
options.DisplayWin = 1;
options.GnFigs = 0;


[posterior,out] = VBA_NLStateSpaceModel(...
    y(:)',u(:)',f_fname,g_fname,dim,options);


% mu = zeros(dim.n_phi,p);
% va = zeros(dim.n_phi,p);
% hf = figure;
% ha = gca(hf);
% set(ha,'nextplot','add')
% for t=1:p
%     [posterior,out] = VBA_NLStateSpaceModel(...
%         y(1:t),u(1:t),[],g_fname,dim,options);
%     mu(:,t) = posterior.muPhi;
%     va(:,t) = diag(posterior.SigmaPhi);
%     if t > 1
%         cla(ha)
%         [haf,hf] = plotUncertainTimeSeries(mu(:,1:t),sqrt(va(:,1:t)),1:t,ha,1:2);
%     end
% end


%---- Display results ----%
displayResults(posterior,out,y,x,x0,theta,[],[],[]);
