% This demo reproduces the variant of the "bookbags and poker chips" task
% of [D'Acremont et al., submitted].
% In brief, it looks at martingale properties of updates of bayesian model
% comparison, where two bernouilli models are compared on the basis of
% random samples of coloured marbles.
%
% At each trial t, a variable u(t) is sampled from a bernouilli
% distribution with probability q0. This variable selects one of two urns,
% from which another variable y(t) is sampled. NB: the probability r0 that
% y(t) is orange depends on u(t), i.e.: r0=P(y(t)=1|u(t)). More precisely,
% r0=0.7 if u(t)=1 and r0=0.3 if u(t)=0. 
% NB: A canonical simulation runs for T trials. There are N repetitions of
% such canonical Monte-Carlo simulations. At each series of trials, q0 was
% drawn from a uniform distribution on the interval [a,b]. 
% A bayesian observer updates its posterior distribution over q0 and u(t). 
% At each trial t, this bayesian observer can predict the next outcome
% y(t+1) according to:
%       p(y(t+1)|y(1,...,t)) = E[P(y(t+1)=1|u(t+1))|y(1,...,t)]
% Let R(t) be the right-hand side of the above equation:
%       R(t) = 0.7*E[q0|y(1,...,t)] + 0.3*(1-E[q0|y(1,...,t)])
% We then look at the following Monte-Carlo expectations:
%  - E[R(t+1),R(t)]
%  - E[(R(t+1)-R(t))*(0.5-R(t))]

clear all
close all


N = 128; % number of Monte-Carlo repetitions
T = 9;   % number of trials or time samples
a = 0; % lower bound of the q0 interval
b = 1; % upper bound of the q0 interval
r0 = [0.7;0.3]; % bernouilli parameters of each urn
v = 1e-2; % concentration parameter of the beta distribution

options.DisplayWin = 0;
c = (a+b)/2;
options.priors.a = 1e-2*[c*v;(1-c)*v];
y = zeros(N,T);
r = zeros(N,T+1);
dr = zeros(N,T+1);
r(:,1) = rand(1,N);

for i=1:N
    q0 = a+rand*(b-a);
    u = sampleFromArbitraryP([q0,1-q0],[1,0],1);
    L = zeros(2,T);
    for t=1:T
        if u
            y = sampleFromArbitraryP([r0(1),1-r0(1)],[1,0],1);
        else
            y = sampleFromArbitraryP([r0(2),1-r0(2)],[1,0],1);
        end
        L(1,t) = y*log(r0(1)) + (1-y).*log(1-r0(1));
        L(2,t) = y*log(r0(2)) + (1-y).*log(1-r0(2));
        [posterior,out] = VBA_groupBMC(L(:,1:t),options);
        r(i,t+1) = r0(1)*out.Ef(1) + r0(2)*(1-out.Ef(1));
        dr(i,t) = r(i,t+1) - r(i,t);
    end
end


mc = zeros(1,T);
sc = zeros(1,T);
ec = zeros(1,T);
m = zeros(1,T);
s = zeros(1,T);
for t=1:T
    mc(t) = mean(dr(:,t).*(0.5-r(:,t)));
    sc(t) = std(dr(:,t).*(0.5-r(:,t)))./sqrt(N);
    m(t) = mean(dr(:,t));
    s(t) = std(dr(:,t))./sqrt(N);
end

hf = figure('color',[1 1 1]);
ha = subplot(2,1,1,'parent',hf,'nextplot','add');
errorbar(ha,m,1.96*s)
plot(ha,[1 T],[0 0],'r')
plot(ha,m,'.')
title(ha,'E[r(t+1)-r(t)]')
legend({'Monte-Carlo simulations','martingale theory'})
ha = subplot(2,1,2,'parent',hf,'nextplot','add');
errorbar(ha,mc,1.96*sc)
plot(ha,[1 T],[0 0],'r')
plot(ha,mc,'.')
legend({'Monte-Carlo simulations','martingale theory'})
title(ha,'Cov[r(t+1)-r(t),0.5-r(t)]')





    
    
    
    