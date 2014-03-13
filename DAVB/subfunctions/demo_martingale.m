% This demo looks at martingale properties of bayesian updates on the
% unique parameter of a binomial likelihood.
% At each trial t, data y(t) is sampled from a bernouilli distribution with
% probability f=1/2. A bayesian observer updates its (beta) posterior
% distribution p(f|y(1,...,t)) over f. A canonical simulation runs for T
% trials. There are N repetitions of such canonical Monte-Carlo
% simulations.
% At each trial t, this bayesian observer can predict the next outcome
% y(t+1) according to:
%       p(y(t+1)|y(1,...,t)) = E[f|y(1,...,t)]
% Let r(t) be the right-hand side of the above equation:
%       r(t) = ( sum(y(1,...,t)) +1)./(t+1)
% We then look at the following Monte-Carlo expectations:
%  - E[r(t+1),r(t)]
%  - E[(r(t+1)-r(t))*(0.5-r(t))]
% Note: here, we have assumed that the Bayesian observer starts with flat
% priors on f, i.e. p(f)=1. This corresponds to setting the parameters of
% the corresponding beta probability density function to 1 (i.e. 1 "dummy"
% count for y=1 and 1 for y=0).

clear all
close all


N = 2^10; % number of Monte-Carlo repetitions
T = 32;   % number of trials or time samples
f = 0.85; % hidden bernouilli frequency

%-- perform Monte-Carlo simulations
y = zeros(N,T);
r = zeros(N,T+1);
dr = zeros(N,T+1);
r(:,1) = 0.5; % prior mean (under flat prior)
for i=1:N
    for t=1:T
        % First, sample from target Bernouilli distribution
        y(i,t) = randn > erfcinv(2*f)*sqrt(2);
        % Then, let bayesian observer assimilate data
        r(i,t+1) = (sum(y(i,1:t))+1)./(t+2);
        % Lastly, store change in update
        dr(i,t) = r(i,t+1) - r(i,t);
    end
end

%-- derive Monte-Carlo averages
mc = zeros(1,T);
sc = zeros(1,T);
m = zeros(1,T);
s = zeros(1,T);
for t=1:T
    mc(t) = mean(dr(:,t).*(0.5-r(:,t)));
    sc(t) = std(dr(:,t).*(0.5-r(:,t)))./sqrt(N);
    m(t) = mean(dr(:,t));
    s(t) = std(dr(:,t))./sqrt(N);
end

hf = figure('color',[1 1 1]);
ha = subplot(3,1,1,'parent',hf,'nextplot','add','xlim',[0,T+0.5]);
mr = mean(r,1);
sr = std(r,[],1)./sqrt(N);
my = mean(y,1);
sy = std(y,[],1)./sqrt(N);
errorbar(ha,[0:T],mr,1.96*sr)
errorbar(ha,my,1.96*sy,'color',[1 0 0])
plot(ha,[0:T],mr,'.')
plot(ha,my,'r.')
title(ha,'E[r(t)]')
legend({'bernouilli parameter estimate','data samples'})
ha = subplot(3,1,2,'parent',hf,'nextplot','add','xlim',[0,T+0.5]);
errorbar(ha,m,1.96*s)
plot(ha,[1 T],[0 0],'r')
plot(ha,m,'.')
title(ha,'E[r(t+1)-r(t)]')
legend({'Monte-Carlo simulations','martingale theory'})
ha = subplot(3,1,3,'parent',hf,'nextplot','add','xlim',[0,T+0.5]);
errorbar(ha,mc,1.96*sc)
plot(ha,[1 T],[0 0],'r')
plot(ha,mc,'.')
legend({'Monte-Carlo simulations','martingale theory'})
title(ha,'Cov[r(t+1)-r(t),0.5-r(t)]')





    
    
    
    