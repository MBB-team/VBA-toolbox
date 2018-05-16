function [posterior, out] = demo_dynamicalSystem ()
% // VBA toolbox //////////////////////////////////////////////////////////
%
% [posterior, out] = demo_dynamicalSystem()
% Demo of dynamical system simulation and inversion
%
%
% Background:
% ~~~~~~~~~~~
%
% /////////////////////////////////////////////////////////////////////////

%% Define the model
% =========================================================================

% Description of the dynamics
% -------------------------------------------------------------------------
% evolution is a simple parametric convolution of the input
f_fname = @f_alpha;
% states are directly observable (y = x + noise)
g_fname = @g_Id;

% evolution function has to be discretized over time:
% - define size of a time step
dt = 1e-1; % in sec
% - store in the structure that will be passed to the evolution function
options.inF.dt = dt;

% Dimensions of the model
% -------------------------------------------------------------------------
% evolution parametesr
dim.n_theta = 2;
% observation parameters
dim.n_phi = 0;
% number of states
dim.n = 1;

%% Simulate data
% =========================================================================

% Define design, here inputs to be convolved
% -------------------------------------------------------------------------
% number of observations
n_t = 1e4 / dt; 
% baseline
u = zeros (1, n_t);
% random pulses
nPulses = 16;
u(randperm (n_t, nPulses)) = 1;

% Parameters of the model to be simulated
% -------------------------------------------------------------------------
% evolution parameters
theta = [1 / dt; 1];
% observation parameters
phi = [];
% initial state
x0 = [0;0];

% state noise precision
alpha = Inf; % deterministic  
% observation precision
sigma = Inf; % no noise

% Simulate
% -------------------------------------------------------------------------
[y, x] = simulateNLSS (n_t, f_fname, g_fname, theta, phi, u, alpha, sigma, options, x0);

% Display
% -------------------------------------------------------------------------
hf = figure('color',[1 1 1]);

timeline = (0 : n_t - 1) * dt;

% plot inputs
ha(1) = subplot (2, 1, 1, 'parent', hf);
plot (ha(1), timeline, u);
xlabel (ha(1), 'time (sec)');
ylabel (ha(1), 'u(t)');
title (ha(1), 'input u');

% plot state trajectory
ha(2) = subplot (2, 1, 2, 'parent', hf);
plot (ha(2), timeline, x);
xlabel (ha(2), 'time (sec)');
ylabel (ha(2), 'x(t)');
title (ha(2), 'output x');

% make pretty
set(ha, 'box', 'off', 'ygrid', 'on', 'ylim', [- 0.2, 1.2]);

%% Perform model estimation, the BAD way
% =========================================================================

% show likelihood
yg = [-0.2:1e-2:1.2];
LL = zeros(length(yg),n_t);
for t=1:n_t
    LL(:,t) = exp(-32*(yg(:)-y(1,t)).^2);
end
hf = figure('color',[1 1 1]);
ha = subplot(1,2,1,'parent',hf,'nextplot','add');
imagesc(LL,'parent',ha)
set(ha(1),'box','off','xlim',[0,n_t-1]);%,'ylim',[-0.2,1.2]);
axis(ha(1),'off')
xlabel(ha(1),'time')
ylabel(ha(1),'p(y|P,m)')
title(ha(1),['likelihood'])

gp1 = theta(1)+[-4:4];
gp2 = theta(2)+[-2:0.5:2];
n1 = length(gp1);
n2 = length(gp2);
for i=1:n1
    i
    for j=1:n2
        theta = [gp1(i);gp2(j)];
        [gy] = simulateNLSS(n_t,f_fname,g_fname,theta,[],u,alpha,sigma,options,x0);
        LLp(i,j) = sum((gy(1,:)-y(1,:)).^2);
    end
end
LLp(find(isnan(LLp)==1)) = Inf;
LLp(LLp>1e3) = 1e3;
ha(2) = subplot(1,2,2,'parent',hf,'nextplot','add');
imagesc(exp(-1e-2*LLp),'parent',ha(2))
set(ha(2),'box','off');
axis(ha(2),'off')
xlabel(ha(2),'time')
ylabel(ha(2),'p(y|P,m)')
title(ha(2),['likelihood'])
      

options.priors.muTheta = .5*ones(2,1);
dim.n = 2;
[posterior,out] = VBA_NLStateSpaceModel(y,u,f_fname,g_fname,dim,options);



