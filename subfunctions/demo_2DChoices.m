% This demo performs a graphical representation of 2D choice data.
% This script simulates an agent making choices between two options, which
% differ in terms of the expected reward and temporal horizons, under a
% hyperbolic temporal discounting model. It then estimates the 2D utility
% profile, using 2D-Fourier basis functions. This semi-parametric estimate
% is then compared to the simulated profile.


clear all
close all

% initialize simulations
ntrials = 1024;
g_fname = @g_discounting;
phi = [-1;0];
in.ind.logk = 1; % index of discount factor (phi)
in.ind.logb = 2; % index of behavioural temperature (phi)
in.ind.t = [1;3]; % index of temporal horizon (u)
in.ind.R = [2;4]; % index of expected reward (u)
in.model = 'hyperbolic';

% set choice alternatives
R = 10*rand(2,ntrials);
t = 10*rand(2,ntrials);
u = zeros(4,ntrials);
u(in.ind.t,:) = t;
u(in.ind.R,:) = R;

% simulate choices under parametric model
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
dr = R(1,:) - R(2,:); % 1st-alternative reward minus 2nd-alternative reward
dt = t(1,:) - t(2,:); % 1st-alternative delay minus 2nd-alternative delay

br = [quantile(dr,.05),quantile(dr,.95)];
bt = [quantile(dt,.05),quantile(dt,.95)];
xr = linspace(br(1),br(2),32);
xt = linspace(bt(1),bt(2),32);
ne  = hist2(dr,dt,xr,xt); % # choice alternatives that fall in that particular domain

hf = figure('color',[1 1 1]);
ha = subplot(2,2,1,'parent',hf,'nextplot','add');
xlabel(ha,'dr')
ylabel(ha,'dt')
ha(2) = subplot(2,2,2,'parent',hf,'nextplot','add');
xlabel(ha(2),'dr')
ylabel(ha(2),'dt')
counts = zeros(size(ne)); % # choices in favour of the 1st-alternative
col = {'r','g'};
for i=1:ntrials
    plot(ha(1),dr(i),dt(i),[col{y(i)+1},'+'])
    counts = counts + (2*y(i)-1)*hist2(dr(i),dt(i),xr,xt);
end
f = counts./(ne+1); % frequency of the 1st-alternative (for each dr and dt)
plot(ha(1),br,[bt(1),bt(1)],'k');
plot(ha(1),br,[bt(2),bt(2)],'k');
plot(ha(1),[br(1),br(1)],bt,'k');
plot(ha(1),[br(2),br(2)],bt,'k');
axis(ha(1),'equal')
axis(ha(1),'tight')
imagesc(f,'parent',ha(2))
axis(ha(2),'equal')
axis(ha(2),'tight')

% invert parametric model
options.priors.SigmaPhi = 1e0*eye(dim.n_phi);
options.DisplayWin = 0;
options.verbose = 1;
[posterior,out] = VBA_NLStateSpaceModel(y,u,[],g_fname,dim,options);
VBA_ReDisplay(posterior,out,1)


% evaluate utility profile on a grid
g1 = linspace(min(vec(R)),max(vec(R)),1e2);
g2 = linspace(min(vec(t)),max(vec(t)),5e1);
N = 8;
X = create2dbf(g1,g2,N);
inb.ind.x = in.ind.t;
inb.ind.y = in.ind.R;
inb.gx = g1;
inb.gy = g2;
inb.bf = X;
g_fname = @g_2AFC_basis;
dim = [];
dim.n_phi = size(X,3);
dim.n = 0;
dim.n_theta = 0;
options.priors.SigmaPhi = 1e0*eye(dim.n_phi);
options.checkGrads = 0;
options.inG = inb;
options.binomial = 1;
options.DisplayWin = 0;
options.verbose = 1;
Eu = zeros(size(X,1),size(X,2));
[p0,o0] = VBA_NLStateSpaceModel(y,u,[],g_fname,dim,options);
VBA_ReDisplay(p0,o0,1)
for ii=1:size(X,3)
    Eu = Eu + X(:,:,ii)*p0.muPhi(ii);
end
nsamples = zeros(size(X,1),size(X,2));
for t=1:ntrials
    u1 = u(inb.ind.x,t);
    u2 = u(inb.ind.y,t);
    for i=1:2
        [tmp,ix] = min((inb.gx-u1(i)).^2);
        [tmp,iy] = min((inb.gy-u2(i)).^2);
        nsamples(ix,iy) = nsamples(ix,iy) + 1;
    end
end

% to be compared with simulated utility function
v = zeros(size(X,1),size(X,2));
for i=1:length(inb.gx)
    for j=1:length(inb.gy)
        uxy = [inb.gx(i);inb.gy(j)];
        [v(i,j)] = v_discounting([],phi,uxy,in);
    end
end
        
        

hf = figure('color',[1 1 1]);
ha = subplot(2,2,1,'parent',hf);
imagesc(Eu,'parent',ha)
title(ha,'estimated utility')
xlabel(ha,'delay')
ylabel(ha,'reward')
ha = subplot(2,2,2,'parent',hf);
imagesc(nsamples,'parent',ha)
title(ha,'# samples')
xlabel(ha,'delay')
ylabel(ha,'reward')
ha = subplot(2,2,3,'parent',hf);
imagesc(v,'parent',ha)
title(ha,'simulated utility')
xlabel(ha,'delay')
ylabel(ha,'reward')
ha = subplot(2,2,4,'parent',hf);
plot(Eu(:),v(:),'.','parent',ha)
title(ha,'comp')
xlabel(ha,'estimated utility')
ylabel(ha,'simulated utility')
getSubplots


