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
EP = sigmoid(m0./(1+a*sqrt(V0)));
VP = EP.*(1-EP).*(1-1./(1+a*sqrt(V0)));

ha = subplot(2,1,2,'parent',hf);
plot(ha,EP','marker','.')
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
