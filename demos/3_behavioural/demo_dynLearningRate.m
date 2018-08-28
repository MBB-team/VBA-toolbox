% demo operant learning with dynamic learning rate
% This demo simulates a "volatile learner" experiencing feedback to her
% chosen actions. The pseudo learning rate then changes as a function of
% (estimated) environmental volatility. We then fit the behavioural
% response with a stochastic variant of the Q-learning model, which allows
% estimating the time series of learning rates. We then perform a Volterra
% decomposition of hidden states (which include the learning rate).
% Finally, we fit an "augmented" Q-learning model, which predicts an
% transient acceleration of learning rate following changesin the winning
% action.


close all
clear variables

% simulate VB volatile learner in a 2-armed bandit task
fb.inH.u0 = repmat([ones(1,50),zeros(1,50)],1,4); % feedbacks with reversals
nt = size(fb.inH.u0,2)+1;
fb.h_fname = @h_truefalse;
fb.indy = 1;
fb.indfb = 2;
theta = [0;-2;0];
phi = [3;0]; % inverse temperature & bias
x0 = repmat([0;0;0;0;0],2,1);
inF.lev2 = 1; % 3rd level (volatility learning)
inF.kaub = 1.4;
inF.thub = 1;
inF.rf = 1;
inG.respmod = 'taylor';
options.sources.type = 1 ; % binomial observation;
options.inF = inF;
options.inG = inG;
options.skipf = zeros(1,nt);
options.skipf(1) = 1;
[y,x,x0,eta,e,u] = VBA_simulate (nt,@f_OpLearn,@g_VBvolatile0,theta,phi,zeros(2,nt),Inf,Inf,options,x0,fb);
% plot simulated behaviour
hf = figure('color',[1 1 1 ],'name','simulated choices');
ha = axes('parent',hf);
plot(ha,y-e,'r')
hold(ha,'on')
plot(ha,y,'kx')
legend(ha,{'p(y=1|theta,phi,m)','binomial data samples'})
VBA_getSubplots ();
dummy.options = options;
[ha,hf] = unwrapVBvolatileOTO(struct('muX',x,'muTheta',theta),dummy);
set(hf,'name','simulated volatile VB learner')


% dummy VB inversion (with ideal priors) of volatile learner
d00 = struct('n',2*5,'n_theta',3,'n_phi',2);
priors = [];
priors.muPhi = phi;
priors.muTheta = theta;
priors.muX0 = x0;
priors.SigmaPhi = 0*eye(d00.n_phi);
priors.SigmaTheta = 0*eye(d00.n_theta);
priors.SigmaX0 = 0*eye(d00.n);
priors.a_alpha = Inf;
priors.b_alpha = 0;
opt00 = options;
opt00.priors = priors;
[p00,o00] = VBA_NLStateSpaceModel(y,u,@f_OpLearn,@g_VBvolatile0,d00,opt00);



%----
% Now set up Q-learner model
u(fb.indfb,u(fb.indfb,:)==0)=-1; % mean-centre feedback for value learning
d0 = struct('n',3,'n_theta',0,'n_phi',1);
priors = [];
priors.a_alpha = 1;
priors.b_alpha = 1;
tmp = 1e2*eye(3);
tmp(3,3) = 1e0;
for t=1:nt
    priors.iQx{t} = tmp;
end
opt0 = [];
opt0.backwardLag = 32;
opt0.priors = priors;
opt0.sources.type = 1;
opt0.verbose = 1;
opt0.MaxIter = 3;
opt0.kernelSize = 32;
opt0.detrendU = 4;
[p0,o0] = VBA_NLStateSpaceModel(y,u,@f_Qlearn_dynLR,@g_softmax,d0,opt0);


% check relation between identified learning rate and volatility of VB-learner:
it = 1:400;
X = [VBA_vec(p0.muX(3,:)),ones(nt,1)];
Y = VBA_vec(sum(x([4,9],:),1));
[pv,stat,df,all] = GLM_contrast(X,Y,[1;0],'F',1,{'learning rate','CST'},{'volatility'});


% perform posterior Volterra decomposition
uu = zeros(3,nt);
uu(1,:) = 2*u(fb.indy,:)-1; % previous own action
for t=1:size(u,2)
    % uu(2,:) = winning action
    if u(fb.indfb,t)==1
        uu(2,t) = uu(1,t);
    else
        uu(2,t) = -uu(1,t);
    end
    % uu(3,:) = winning action stability
    try
        if ~isequal(uu(2,t),uu(2,t-1))
            uu(3,t) = 1;
        else
            uu(3,t) = 0;
        end
    end
end
o0_v = o0;
o0_v.u = uu;
[kernels] = VBA_getVolterraKernels(p0,o0_v);
hf = figure('color',[1 1 1],'name','Volterra decomposition');
mk = squeeze(kernels.x.m(3,:,:))';
vk = squeeze(kernels.x.v(3,:,:))';
ha = subplot(2,1,1,'parent',hf);
plotUncertainTimeSeries(mk,vk,[],ha)
legend(ha,{'agent''s chosen action','winning action','winning action stability'})
xlabel(ha,'lag')
ylabel(ha,'Volterra weight')
title(ha,'stochastic learning rate')



% Now fit agent's choice with augmented Q-learning model
u(3,:) = uu(3,:); % add in winning action stability for learning rate evolution
u(3,u(3,:)==0)=-1; % mean-centre for arbitrary steady-state of learning rate
d1 = struct('n',4,'n_theta',2,'n_phi',1);
priors = [];
priors.a_alpha = Inf;
priors.b_alpha = 0;
opt1.priors = priors;
opt1.sources.type = 1;
opt1.figName = 'augmented Q-learning model';
[p1,o1] = VBA_NLStateSpaceModel(y,u,@f_Qlearn_gammaLR,@g_softmax,d1,opt1);


% extract and plot Volterra kernels
o1_v = o1;
o1_v.u = uu;
[kernels] = VBA_getVolterraKernels(p1,o1_v);
mk = squeeze(kernels.x.m(3,:,:))';
vk = squeeze(kernels.x.v(3,:,:))';
ha = subplot(2,1,2,'parent',hf);
plotUncertainTimeSeries(mk,vk,[],ha)
legend(ha,{'agent''s chosen action','winning action','winning action stability'})
xlabel(ha,'lag')
ylabel(ha,'Volterra weight')
title(ha,'augmented learning rate')


