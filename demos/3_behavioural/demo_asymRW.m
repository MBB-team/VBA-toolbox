% Demo: Comparison of asymmetric reinforcement learning models.



close all
clear all

% Basic settings
n_t = 128; % number of trials
f_fname = @f_rwl2; % modified Rescorla-Wagner learning rule
g_fname = @g_goNogo; % -> 'go' choice probability

% simu parameters
feedbacks = [-1; 0; 1]; % feedbacks: negative, neutral, positive
contingencies = [2; 1; 2] / 5; %probability of feedbacks
truemodel = 1; % index of true model (which generates the data)
theta = [2;0;1;0];
phi = [-1]; % value of the 'no-go' option
sigma = 1e1; % precision on value signal
x0 = 0; % initial conditions


% model 1: gain/loss asymmetry in experienced reward
d{1}.n = 1; % hidden states (=value of the 'go' options)
d{1}.n_phi = 1; % observation parameters (bias, temperature)
d{1}.n_theta = 4; % evolution parameters (depends on model)
opt{1}.inF.model = 'utility';
opt{1}.inF.indR = [1:3]; % feedback weigths
opt{1}.inF.indAlpha = 4; % learning rate
opt{1}.inF.indu = 1; % feedback
opt{1}.inF.inda = 2; % last choice
opt{1}.priors.a_alpha = Inf; % deterministic system
opt{1}.priors.b_alpha = 0; % [id]
opt{1}.sources.type = 1; % binary (go/nogo) choices
opt{1}.skipf = zeros(1,n_t);
opt{1}.skipf(1) = 1; % apply identity mapping from x0 to x1.


% model 2: gain/loss asymmetry in learning rate
d{2} = d{1};
opt{2} = opt{1};
opt{2}.inF.model = 'learning';
opt{2}.inF.indR = 1; % feedback weight
opt{2}.inF.indAlpha = [2:4]; % learning rates

% model 3: both asymmetries
d{3} = d{1};
d{3}.n_theta = 6;
opt{3} = opt{1};
opt{3}.inF.model = 'both';
opt{3}.inF.indR = [1:3]; % feedback weigths
opt{3}.inF.indAlpha = [4:6]; % learning rates

% model 4: no asymmetry
d{4} = d{1};
d{4}.n_theta = 2;
opt{4} = opt{1};
opt{4}.inF.model = 'none';
opt{4}.inF.indR = 1; % feedback weigths
opt{4}.inF.indAlpha = 2; % learning rates

% sample random feedback sequence
u = VBA_random ('Arbitrary', contingencies, feedbacks, 1, n_t);

% simulate agent's response given predefined feedback sequence
fb.h_fname = @h_goNogo;
fb.indy = 2;
fb.indfb = 1;
fb.inH.u = u;
[y_choice,x,x0,eta,ee,uu] = VBA_simulate (n_t,f_fname,g_fname,theta,phi,[0;0],Inf,[],opt{truemodel},x0,fb);
y_value = x + randn(size(x))./sqrt(sigma);

U = zeros(3,n_t);
for t=1:n_t
    if y_choice(t)
        iu = u(1,t)+2; % 1 if negative, 2 if neutral, 3 if positive
        U(iu,t) = 1;
    else
        U(2,t) = .2;
    end
end
hf = figure('color',[1 1 1]);
ha = subplot(3,2,1,'box','off','parent',hf);
imagesc(U,'parent',ha)
title(ha,'feedbacks')
xlabel(ha,'time (trials)')
hold(ha,'on')
plot(get(ha,'xlim'),[1.5 1.5],'k')
plot(get(ha,'xlim'),[2.5 2.5],'k')
set(ha,'ytick',[1,2,3],'yticklabel',{'negative','neutral','positive'})
colormap(flipud(bone))

ha = subplot(3,2,2,'box','off','parent',hf);
hp = plot(ha,x,'r');
set(hp,'LineWidth',2)
hold(ha,'on')
plot(ha,y_choice,'k.')
plot(ha,y_value,'k')
legend(ha,{'true value','choice','measured value'})
set(ha,'xlim',[1,n_t])
box(ha,'off')
title(ha,'agent''s response')
xlabel(ha,'time (trials)')


% invert all models
str = cell(4,1);
for i=1:4
    opt{i}.DisplayWin = 1;
    % invert models on choice data
    opt{i}.figName = [opt{i}.inF.model, ' (choice data)'];
    opt{i}.sources.type = 1;
    [p{i,1},o{i,1}] = VBA_NLStateSpaceModel(y_choice,uu,f_fname,g_fname,d{i},opt{i});
    F(i,1) = o{i,1}.F;
    % invert models on value data
    opt{i}.figName = [opt{i}.inF.model, ' (value data)'];
    opt{i}.sources.type = 0;
    [p{i,2},o{i,2}] = VBA_NLStateSpaceModel(y_value,uu,f_fname,@g_Id,d{i},opt{i});
    F(i,2) = o{i,2}.F;
    % store for future display
    str{i} = opt{i}.inF.model;
end

pm(:,1) = VBA_softmax(F(:,1));
pm(:,2) = VBA_softmax(F(:,2));
pf1(:,1) = [pm(1,1)+pm(3,1);pm(2,1)+pm(4,1)];
pf1(:,2) = [pm(1,2)+pm(3,2);pm(2,2)+pm(4,2)];
pf2(:,1) = [pm(2,1)+pm(3,1);pm(1,1)+pm(4,1)];
pf2(:,2) = [pm(2,2)+pm(3,2);pm(1,2)+pm(4,2)];

ha = subplot(3,2,3,'box','off','parent',hf);
hb = bar(ha,pm(:,1),'r');
set(hb,'facecolor',0.8*ones(1,3))
hold(ha,'on')
plot(ha,[0 5],[0.95 0.95],'r--')
set(ha,'xlim',[0,5],'xtick',[1:4],'xticklabels',str)
box(ha,'off')
title(ha,'asymmetry in choice data')
xlabel(ha,'models')
ylabel(ha,'P(m|y,u)')

ha = subplot(3,2,4,'box','off','parent',hf);
hb = bar(ha,pm(:,2),'r');
set(hb,'facecolor',0.8*ones(1,3))
hold(ha,'on')
plot(ha,[0 5],[0.95 0.95],'r--')
set(ha,'xlim',[0,5],'xtick',[1:4],'xticklabels',str)
box(ha,'off')
title(ha,'asymmetry in value data')
xlabel(ha,'models')
ylabel(ha,'P(m|y,u)')


ha = subplot(3,4,9,'box','off','parent',hf);
hb = bar(ha,pf1(:,1),'r');
set(hb,'facecolor',0.8*ones(1,3))
hold(ha,'on')
plot(ha,[0 5],[0.95 0.95],'r--')
set(ha,'xlim',[0,3],'xtick',[1:2],'xticklabels',{'yes','no'})
box(ha,'off')
title(ha,'choice data')
xlabel(ha,'asymmetric utility')
ylabel(ha,'P(f|y,u)')

ha = subplot(3,4,11,'box','off','parent',hf);
hb = bar(ha,pf1(:,2),'r');
set(hb,'facecolor',0.8*ones(1,3))
hold(ha,'on')
plot(ha,[0 5],[0.95 0.95],'r--')
set(ha,'xlim',[0,3],'xtick',[1:2],'xticklabels',{'yes','no'})
box(ha,'off')
title(ha,'value data')
xlabel(ha,'asymmetric utility')
ylabel(ha,'P(f|y,u)')


ha = subplot(3,4,10,'box','off','parent',hf);
hb = bar(ha,pf2(:,1),'r');
set(hb,'facecolor',0.8*ones(1,3))
hold(ha,'on')
plot(ha,[0 5],[0.95 0.95],'r--')
set(ha,'xlim',[0,3],'xtick',[1:2],'xticklabels',{'yes','no'})
box(ha,'off')
title(ha,'choice data')
xlabel(ha,'asymmetric learning')
ylabel(ha,'P(f|y,u)')

ha = subplot(3,4,12,'box','off','parent',hf);
hb = bar(ha,pf2(:,2),'r');
set(hb,'facecolor',0.8*ones(1,3))
hold(ha,'on')
plot(ha,[0 5],[0.95 0.95],'r--')
set(ha,'xlim',[0,3],'xtick',[1:2],'xticklabels',{'yes','no'})
box(ha,'off')
title(ha,'value data')
xlabel(ha,'asymmetric learning')
ylabel(ha,'P(f|y,u)')

drawnow

% for i=1:4
%     VBA_ReDisplay(p{i},o{i},1);
% end



