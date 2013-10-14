function demo_Qlearning

% demo Q-learning
% This demo simulates a number of sequences of choices of a Q-learning
% agent, and then estimates its expected Volterra kernels, given a basis
% function set that captures the susceptibility to the agent's own
% feedback, as well as the susceptibility to the agent's opponent's
% feedback.

% close all
% clear variables
% clc

Nmcmc = 64;

% simulation parameters
theta = sigm(0.5,struct('INV',1)); % learning rate = 0.75
phi = log(1); % inverse temperature = 2

f_fname = @f_Qlearn2;
g_fname = @g_softmax;
h_fname = @h_randOutcome;%@h_truefalse;

fb.inH.er = 1;
fb.inH.vr = 0;
fb.h_fname = h_fname;
fb.indy = 1;
fb.indfb = 2;

dim = struct('n',2,'n_theta',1,'n_phi',1);

priors.muPhi = zeros(dim.n_phi,1);
priors.muTheta = zeros(dim.n_theta,1);
priors.muX0 = zeros(2,1);
priors.SigmaPhi = 1e1*eye(dim.n_phi);
priors.SigmaTheta = 1e1*eye(dim.n_theta);
priors.SigmaX0 = 1e1*eye(dim.n);
priors.a_alpha = Inf;
priors.b_alpha = 0;
options.priors = priors;
options.binomial = 1;
options.verbose = 0;

tau = 8;

nt = 100;
ind = 2:6*nt;

for i=1:Nmcmc
    
    
    
    % allocate feedback struture for simulations
    u0 = [randn(1,nt)>-0.05]; % possible feedbacks
    fb.inH.u0 = [u0,~u0,u0,~u0,u0,~u0]; % with reversals
    
    % choose dummy initial conditions
    x0 = zeros(2,1);
    u = zeros(2,size(fb.inH.u0,2)+1);
    n_t = length(u); % number of trials
    options.skipf = zeros(1,n_t);
    options.skipf(1) = 1; % apply identity mapping from x0 to x1.
    [y,x,x0,eta,e,u] = simulateNLSS_fb(n_t,f_fname,g_fname,theta,phi,u,Inf,Inf,options,x0,fb);
    
    %     figure
    %     plot(y-e,'r')
    %     hold on
    %     plot(y,'kx')
    %     legend({'p(y=1|theta,phi,m)','binomial data samples'})
    %
    %     [p0,o0] = VBA_NLStateSpaceModel(y,u,f_fname,g_fname,dim,options);
    %
    %     displayResults(p0,o0,y,x,x0,theta,phi,Inf,Inf)
    %
    %     getSubplots
    
    
    % form basis function set
    u1 = u(1,:); % own action
    u3 = u(2,:); % feedback
    u2 = zeros(size(u1)); % opponent's action
    u2(u3>0) = u1(u3>0);
    u2(u3<0) = 1-u1(u3<0);
    uu = 2*[u1;u2]-1;
    uu = uu(:,ind);
    y = y(ind);
    
%     corrcoef(uu')
%     figure,plot(uu')
%     [uu] = VBA_orth(uu',1)';
%     corrcoef(uu')
%     figure,plot(uu')
%         pause
    
    
    
    
    % estimate Volterra kernel
    nu = size(uu,1);
    opt.inG.dim.n_t = tau;
    opt.inG.dim.nu = nu;
    [opt.inG.dgdp] = VBA_conv2glm(uu,tau); % shape convolution GLM
    d.n = 0;
    d.n_theta = 0;
    d.n_phi = tau*nu +1;
    opt.priors.muPhi = zeros(d.n_phi,1);
    opt.priors.SigmaPhi = 1e0*eye(d.n_phi);
    opt.binomial = 1;
    opt.DisplayWin = 0;
    opt.verbose = 0;
    [posterior,out] = VBA_NLStateSpaceModel(y',[],[],@g_convSig,d,opt);
    
    [mw(:,:,i),vw(:,:,i)] = extractKernels(posterior,out);
    
%     hf0 = figure('color',[1 1 1 ]);
%     plotVolterra(hf0,mw(:,:,i),vw(:,:,i))
%     pause
    
    i
    
end

hf = figure('color',[1 1 1]);
ha = plotVolterra(hf,mean(mw,3),var(mw,[],3));

hold(ha(1),'on')
plot(ha(1),exp(phi)*sigm(theta)*(1-sigm(theta)).^[0:tau-1],'marker','.')
legend(ha(1),{'estimate','credible interval','theoretical'})


VBA_ReDisplay(posterior,out,1)

% VBA_ReDisplay(p0,o0,1)


