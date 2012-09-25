% This demo looks at martingale properties of RW updates on the unique
% parameter of a binomial likelihood. 
% At each trial t, data y(t) is sampled from a bernouilli distribution with
% probability r=1/2. A RW observer updates its belief about r. A canonical
% simulation runs for T trials. There are N repetitions of such canonical
% Monte-Carlo simulations.
% At each trial t, this RW observer can predict the next outcome
% y(t+1) according to:
%       p(y(t+1)|y(1,...,t)) = r_RW(t)
% The RW observer's learning rule is as follows:
%       r_RW(t+1) = r_RW(t) + alpha*(y(t+1)-r_RW(t))
% where alpha is the learning rate
% We then look at the following Monte-Carlo expectations:
%  - E[R(t+1),R(t)]
%  - E[(R(t+1)-R(t))*(0.5-R(t))]

clear all
close all


N = 2^11; % number of Monte-Carlo repetitions
T = 32;   % number of trials or time samples
alpha = 0.75;

y = zeros(N,T);
r = zeros(N,T+1);
dr = zeros(N,T+1);
r(:,1) = rand(1,N);
for i=1:N
    for t=1:T
        y(i,t) = randn>0.75;
        r(i,t+1) = r(i,t) + alpha*(y(i,t)-r(i,t));
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
    B = [1;alpha]*var(0.5-r(:,t))*[1;alpha]';
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





    
    
    
    