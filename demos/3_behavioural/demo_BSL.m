function demo_BSL
% This script simulates and inverts a Bayesian sequence learner (BSL)

% simulation parameters
K = 1;
options.inF = struct('K',K);
options.inG = struct('K',K);
dim.n = 2^(K+1);
dim.n_theta = 1;
dim.n_phi = 2;



%% simulate sequence of BSL choices
x0 = 1*ones(dim.n,1); % log-odds of P(o=1)
theta = [-2]; % BSL's prior volatility
phi = [-1;0]; % (log-) inverse-temperature, bias
N = 150; % number of trials
p = 0.75;
P = repmat([p,1-p,1-p,p],1,N); % probabilistic repetition of [1 0 0 1]
% P = repmat(p,1,N); % probabilistic repetition of [1 0 0 1]
y = VBA_random ('Bernoulli', P(1 : N));

a = zeros(1,N);
gx = NaN(1,N);
x = zeros(dim.n,N);
x(:,1:K) = repmat(x0,1,K); %initialize hidden states
u(:,1:K) = NaN(K+1,K);
for i=K+2:N
    u(:,i) = flipud(VBA_vec(y(i-K-1:i-1)));
    if K==0 && i==1 % the issue only arises for 0-BSL (degenerated!)
        x(:,i) = f_BSL(x0,theta,u(:,i),options.inF);
    else
        x(:,i) = f_BSL(x(:,i-1),theta,u(:,i),options.inF);
    end
    gx(i) = g_BSL(x(:,i),phi,u(:,i),options.inG);
    a(i) = VBA_random ('Bernoulli', gx(i));
end
hf = figure('color',[1 1 1]);
ha = subplot(2,1,1,'parent',hf);
plot(ha,x')
ha = subplot(2,1,2,'parent',hf);
plot(ha,gx)
hold on, plot([1,N],[0.5,0.5],'color',0.2*[1 1 1])
ic = find(a(1:N)==y(1:N));
plot(ha,ic,a(ic),'g.')
ii = find(a(1:N)~=y(1:N));
plot(ha,ii,a(ii),'r.')
title(ha,['perf=',num2str(100*length(ic)./N,3),'%'])


%% invert BSL model given sequence of agent's choices
options.skipf = zeros(1,N);
options.skipf(1:K+1) = 1;
options.sources = struct('type',1 ,'out', 1); % one binomial observation;
options.SigmaTheta = 1e2*eye(dim.n_theta);
f_fname = @f_BSL;
g_fname = @g_BSL;
[posterior,out] = VBA_NLStateSpaceModel(a,u,f_fname,g_fname,dim,options);

%% Display results
displayResults(posterior,out,a,x,x0,theta,phi,[],[])
hf = unwrapKBSL(posterior.muX,posterior.muPhi,u,options.inG);
end

%% ########################################################################
function hf = unwrapKBSL(x,phi,u,inG)
% display k-BSL's evolving beliefs over trials
% function unwrapKBSL(x,options)
% Bayesian Sequence Learners (BSL) essentially update their posterior
% belief about the transition probabilities of a sequence of outcomes. This
% function displayes BSL's belief as trials unfold. [See f_BSL.m]
% IN:
%   - x: k-BSL's hidden-states
%   - inG: the input structure to g_kBSL.m
% OUT:
%   - hf: display figure handle

K = inG.K; % sequence depth

nt = size(x,2); % nb of trials

hf = figure('color',[1 1 1],'name',[num2str(K),'-BSL learner']);

% derive BSL's prediction about next outcome
gx = 0.5*ones(1,nt);
for t=K+1:nt
    gx(t) = g_BSL(x(:,t),phi,u(:,t),inG);
end
ha = subplot(2,1,1,'parent',hf);
plot(ha,gx,'color','k','marker','.')
xlabel(ha,'time/trials')
ylabel(ha,'P(a=1)')
title(ha,'BSL''s next bet')
box(ha,'off')

% BSL's learned belief about next outcome (given past sequence)
m0 = x(1:2^K,:); % E[log-odds]
V0 = exp(x((2^K)+1:2^(K+1),:)); % V[log-odds]
a = 0.368;
EP = VBA_sigmoid(m0./(1+a*sqrt(V0)));
VP = EP.*(1-EP).*(1-1./(1+a*sqrt(V0)));

ha = subplot(2,1,2,'parent',hf);
plotUncertainTimeSeries(EP,VP,[],ha);
% plot(ha,EP','marker','.')
xlabel(ha,'time/trials')
ylabel(ha,'P(o=1|past o)')
title(ha,'BSL''s conditional belief about P(o|past o)')
box(ha,'off')

if K>0
    leg = cell(2^K,1);
    for k=0:2^K-1
        tmp = dec2bin(k);
        ntmp = length(num2str(tmp));
        leg{k+1} = ['past o = ',cat(2,repmat('0',1,K-ntmp),tmp)];
    end
    legend(ha,leg)
end
end
