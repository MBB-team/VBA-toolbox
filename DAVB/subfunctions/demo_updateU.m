% Demo for utility learning
% This demo simulates how an agent's set of preferences are modified as a
% consequence of the agent's choices, and in the absence of information
% about the consequence of the choices (e.g.: reward). More precisely, we
% reproduce the typical paradigm that has been used to measure the
% phenomenon of 'cognitive dissonance', which consists in three steps:
% 1- people are asked to rate how much they like a set of n items
% 2- people chosee between random paris of items
% 3- people are asked, again, to rate how much they like the items
% The typical result is that the second preference rating of chosen (resp.
% unchosen) items is better (resp. worse), on average, than the first.

clear variables
close all



% %---- First, check how informative are choices
% % NB: this is to evaluate the statistical confound to cognitive dissonance
% % effects
% 
% 
% % 1- assign random utility U0 to the items
% in.n = 16; % # items
% U0 = randn(in.n,1);
% U0 = U0 - mean(U0); % for comparison with estimated utility
% 
% % 2- Measure preference ratings
% % NB: preference ratings are noisy estimates of U0
% sr = 1e0; % noise std in preference ratings
% r1 = U0 + sr*randn(in.n,1);
% r2 = U0 + sr*randn(in.n,1);
% N = 128; % # MCMC simulations (see below)
% 
% % 3- choice task
% % NB: here the choices are based upon the underlying utility function U0
% n_t = 512; % # forced choice trials
% u = zeros(2,n_t);
% beta = 1; % behavioural temperature
% for t=1:n_t
%     % pick a set of 2 items at random
%     ind = randperm(in.n);
%     ind = ind(1:2);
%     % sample choice according to the probability
%     dv = U0(ind(1)) - U0(ind(2)); % relative value of item 1
%     p = sigm(dv/beta); % probability of picking the first item
%     c = sampleFromArbitraryP([p,1-p],[1,0]);
%     % store chosen and unchosen item
%     u(:,t) = c*ind(:) + (1-c)*flipud(ind(:));
% end
% 
% % 4- reconstruct utility profile from choices
% dim.n = 0;
% dim.n_phi = in.n +1;
% dim.n_theta = 0;
% g_fname = @g_u2c;
% priors.SigmaPhi = eye(dim.n_phi);
% priors.SigmaPhi(end,end) = 0; % fix temperature to 1
% options.priors = priors;
% options.inG.temp = in.n+1;
% options.inG.effort = [];
% options.inG.ic = 1;
% options.inG.iu = 2;
% options.binomial = 1;
% [posterior,out] = VBA_NLStateSpaceModel(ones(1,n_t),u,[],g_fname,dim,options);
% 
% % 5- check accuracy of utility profile estimate against simulated truth
% Eu = posterior.muPhi(1:in.n);
% Vu = diag(posterior.SigmaPhi);
% Vu = Vu(1:in.n);
% hf = figure('color',[1 1 1]);
% ha = subplot(3,2,1,'parent',hf,'nextplot','add');
% plotUncertainTimeSeries(Eu,Vu,[],ha)
% set(ha,'xlim',[0,in.n+1])
% plot(ha,U0,'go')
% legend(ha,{'estimated utility (Eu)','credible interval','simulated utility (U0)'})
% xlabel('items')
% ylabel('utility')
% ha = subplot(3,2,2,'parent',hf,'nextplot','add');
% plot(ha,U0,Eu,'b+')
% grid(ha,'on')
% xlabel('simulated utility (U0)')
% ylabel('estimated utility (Eu)')
% 
% 
% ha = subplot(3,2,3,'parent',hf,'nextplot','add');
% plot(ha,U0,r1,'b+')
% grid(ha,'on')
% xlabel('simulated utility (U0)')
% ylabel('first rating (r1)')
% 
% ha = subplot(3,2,4,'parent',hf,'nextplot','add');
% plot(ha,U0,r2,'b+')
% grid(ha,'on')
% xlabel('simulated utility (U0)')
% ylabel('second rating (r2)')
% 
% % 6- try to explain second rating with first rating and utility profile
% % estimated from choices. The idea is that both ratings and the estimated
% % utility profile are noisy proxies for the true utility. Thus, the average
% % of the first rating and and the estimate utility profile should be a good
% % approximation to the second rating...
% R = zeros(4,N);
% for i=1:N
%     r1 = U0 + sr*randn(in.n,1);
%     r2 = U0 + sr*randn(in.n,1);
%     X = r1;
%     [pv,stat,df,all] = GLM_contrast(X,r2,1,'t',0,{'r1'},{'r2'});
%     R(1,i) = all.R2;
%     X = 0.5*(r1+Eu);
%     [pv,stat,df,all] = GLM_contrast(X,r2,1,'t',0,{'r1','choices'},{'r2'});
%     R(2,i) = all.R2;
%     X = Eu;
%     [pv,stat,df,all] = GLM_contrast(X,r2,1,'t',0,{'r1','choices'},{'r2'});
%     R(3,i) = all.R2;
%     X = U0;
%     [pv,stat,df,all] = GLM_contrast(X,r2,1,'t',0,{'r1','choices'},{'r2'});
%     R(4,i) = all.R2;
% end
% 
% ha = subplot(3,2,5,'parent',hf,'nextplot','add');
% mR = mean(R,2);
% vR = var(R,[],2);
% plotUncertainTimeSeries(mR,vR,[],ha)
% set(ha,'xlim',[0 5],'xtick',1:4,'xticklabels',{'r1','Eu','(r1+Eu)/2','U0'})
% xlabel('GLM regressors')
% ylabel('r2 fit quality')
% 




%---- Second, simulate cognitive dissonance
% NB: now, we simulate how choices modify preferences according to the
% 'closed' scenario

% 1- assign flat utility U0 to the items
n = 2; % # items
U0 = 0.*randn(n,1);
U0 = U0 - mean(U0); % for comparison with estimated utility

% 2- Perform first preference rating r1
sr = 1e0; % noise std in preference ratings
r1 = U0 + sr*randn(n,1);

% 3- choice task
% NB: here the choices influence the underlying utility function
n_t = 1024; % # forced choice trials
fb.indfb = 4:5;
fb.indy = 3;
fb.h_fname = @h_whichItem;
fb.inH.ind = zeros(2,n_t);
% pick pairs of items to compare
for t=1:n_t
    ind = randperm(n);
    fb.inH.ind(:,t) = vec(ind(1:2));
end
f_fname = @f_updateU;
g_fname = @g_u2p;
dim.n = n+n^2;
dim.n_phi = 1;
dim.n_theta = 1;
inG.n = n;
inG.temp = 1;
inG.iu = 1:2;
inG.ic = fb.indy;
inG.effort = [];
inF.n = n;
inF.temp = 1;
inF.iu = fb.indfb;
inF.ic = fb.indy;
inF.effort = [];
options.inG = inG;
options.inF = inF;
options.binomial = 1;
options.skipf = zeros(1,n_t);
options.skipf(1) = 1;
u = zeros(5,n_t);
u(inG.iu,:) = fb.inH.ind;
x0(1:n,1) = U0; % initial utility
S0 = eye(n); % initial variance on utility
x0(n+1:n+n^2) = vec(S0); % initial variance
phi = -4; % log temp
theta = 0; % believed log temp
[y,x,x0,eta,e,u] = simulateNLSS_fb(n_t,f_fname,g_fname,0,0,u,Inf,[],options,x0,fb);

mu = x(1:n,:);
figure,plot(mu')


%
%
% % Display results
% displayResults(posterior,out,y,x,x0,theta,phi,alpha,sigma)


