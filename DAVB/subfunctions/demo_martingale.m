% This demo looks at martingale properties of bayesian updates on the
% unique parameter of a binomial likelihood.
% At each trial t, data y(t) is sampled from a bernouilli distribution with
% probability r=1/2. A bayesian observer updates its (beta) posterior
% distribution p(r|y(1,...,t)) over r. A canonical simulation runs for T
% trials. There are N repetitions of such canonical Monte-Carlo
% simulations.
% At each trial t, this bayesian observer can predict the next outcome
% y(t+1) according to:
%       p(y(t+1)|y(1,...,t)) = E[r|y(1,...,t)]
% Let R(t) be the right-hand side of the above equation:
%       R(t) = sum(y(1,...,t))./t
% We then look at the following Monte-Carlo expectations:
%  - E[R(t+1),R(t)]
%  - E[(R(t+1)-R(t))*(0.5-R(t))]

clear all
close all


N = 2^11; % number of Monte-Carlo repetitions
T = 32;   % number of trials or time samples

y = zeros(N,T);
r = zeros(N,T+1);
dr = zeros(N,T+1);
r(:,1) = rand(1,N);
for i=1:N
    for t=1:T
        y(i,t) = randn>0;
        r(i,t+1) = sum(y(i,1:t))./t;
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
    B = [1;1/t]*var(0.5-r(:,t))*[1;1/t]';
    ec(t) = B(2,1);
end

hf = figure('color',[1 1 1]);
ha = subplot(2,1,1,'parent',hf,'nextplot','add');
errorbar(ha,m,1.96*s)
plot(ha,[1 T],[0 0],'r')
plot(ha,m,'.')
title(ha,'E[r(t+1)-r(t)]')
legend({'Monte-Carlo simulations','theory'})
ha = subplot(2,1,2,'parent',hf,'nextplot','add');
errorbar(ha,mc,1.96*sc)
plot(ha,ec,'r')
plot(ha,mc,'.')
legend({'Monte-Carlo simulations','theory'})
title(ha,'Cov[r(t+1)-r(t),0.5-r(t)]')





    
    
    
    