% This demo performs a graphical representation of 2D choice data.
% This script simulates an agent making choices between two options, which
% differ in terms of the expected reward and temporal horizons, under a
% hyperbolic temporal discounting model.  

clear all
close all

% initialize simulations
ntrials = 1024;
g_fname = @g_discounting;
phi = [1;-1];
in.ind.logk = 1; % index of discount factor (phi)
in.ind.logb = 2; % index of behavioural temperature (phi)
in.ind.t = [1;3]; % index of temporal horizon (u)
in.ind.R = [2;4]; % index of expected reward (u)
in.model = 'hyperbolic';

% set choice alternatives
R = exp(randn(2,ntrials));
t = exp(randn(2,ntrials));
u = zeros(4,ntrials);
u(in.ind.t,:) = t;
u(in.ind.R,:) = R;

% simulate choices
dim.n_theta = 0;
dim.n_phi = 2;
dim.n = 0;
dim.p = 1;
dim.n_t = ntrials;
options.inG = in;
options.binomial = 1;
options.dim = dim;
[y,x,x0,eta,e] = simulateNLSS(ntrials,[],g_fname,[],phi,u,[],[],options);


% graphical summary of choice data
dr = R(1,:) - R(2,:);
dt = t(1,:) - t(2,:);

br = [quantile(dr,.05),quantile(dr,.95)];
bt = [quantile(dt,.05),quantile(dt,.95)];
xr = linspace(br(1),br(2),32);
xt = linspace(bt(1),bt(2),32);
ne  = hist2(dr,dt,xr,xt);

hf = figure('color',[1 1 1]);
ha = subplot(2,2,1,'parent',hf,'nextplot','add');
xlabel(ha,'dr')
ylabel(ha,'dt')
ha(2) = subplot(2,2,2,'parent',hf,'nextplot','add');
xlabel(ha(2),'dr')
ylabel(ha(2),'dt')
counts = zeros(size(ne));
col = {'r','g'};
for i=1:ntrials
    plot(ha(1),dr(i),dt(i),[col{y(i)+1},'+'])
    counts = counts + (2*y(i)-1)*hist2(dr(i),dt(i),xr,xt);
end
plot(ha(1),br,[bt(1),bt(1)],'k');
plot(ha(1),br,[bt(2),bt(2)],'k');
plot(ha(1),[br(1),br(1)],bt,'k');
plot(ha(1),[br(2),br(2)],bt,'k');
axis(ha(1),'equal')
axis(ha(1),'tight')

imagesc(counts./(ne+1),'parent',ha(2))
axis(ha(2),'equal')
axis(ha(2),'tight')

% invert model
options.priors.SigmaPhi = 1e0*eye(dim.n_phi);
[posterior,out] = VBA_NLStateSpaceModel(y,u,[],g_fname,dim,options);

out.u = [dr;dt];
[kernels] = VBA_VolterraKernels(posterior,out,16);
out.diagnostics.kernels = kernels;
VBA_ReDisplay(posterior,out,1)





