% demo for split-Laplace (MoG approx to a Gaussian density)

clear variables
close all

% % 0- evaluate the MoG approximation to the N(0,1) density.
% nmog = 2;
% [m,s,p] = getMoG4N01(nmog,1,1);

% 1- evaluate the split-Laplace approximation to the prior predictive
% density.
% NB: use 1D density with sigmoid observation mapping.

g_fname = @g_x2;
dx = 2e-2;
gridy = -20:dx:20;
gridx = -20:dx:20;
ng = length(gridx);
mx = 0;
vx = 1e0;
px = exp(-0.5*(gridx-mx).^2./vx)./sqrt(2*pi*vx);
px = px./sum(px);
vy = 1e0;
py = zeros(size(gridy));
for i=1:ng
    gx = feval(g_fname,[],gridx(i),[],[]);
    tmp = exp(-0.5*(gridy-gx).^2./vy)./sqrt(2*pi*vy);
    py = py + tmp.*px(i);
end
py = py./sum(py);
Ey0 = sum(gridy.*py);
Vy0 = sum((gridy-Ey0).^2.*py);
pyn = -0.5.*(gridy-Ey0).^2./Vy0;
pyn = exp(pyn-max(pyn));
pyn = pyn./sum(pyn);

nmax = 10;
inG.x = 1;
priors.muPhi = 0;
priors.SigmaPhi = vx;
priors.a_sigma = 1;
priors.b_sigma = vy;
options.inG = inG;
options.priors = priors;
options.checkGrads = 0;
dim.n_phi = 1;
dim.n = 0;
dim.n_theta = 0;
dim.p = 1;
dim.n_t = 1;
muy = zeros(nmax,1);
Vy = zeros(nmax,1);
for nmog=1:nmax
    nmog
    [muy(nmog),Vy(nmog),m,V,w] = splitLaplace(...
        [],...
        [],...
        g_fname,...
        dim,...
        options,...
        nmog);
end

% last split-Laplace approximation to the prior predictive density
pymog = -0.5.*(gridy-muy(end)).^2./Vy(end);
pymog = exp(pymog-max(pymog));
pymog = pymog./sum(pymog);

pylap = -0.5.*(gridy-muy(1)).^2./Vy(1);
pylap = exp(pylap-max(pylap));
pylap = pylap./sum(pylap);

hf = figure('color',[1 1 1]);
ha = subplot(2,1,1,'parent',hf);
set(ha,'nextplot','add')
% plot prior predictive density
plot(ha,gridy,py,'r');
% plot moment-matching approximation to the prior predictive density
plot(ha,gridy,pyn,'g');
% plot Laplace approximation to the prior predictive density
plot(ha,gridy,pylap,'k:');
% plot last split-Laplace approximation to the prior predictive density
plot(ha,gridy,pymog,'k--');
set(ha,'xlim',[-10,10])
grid(ha,'on')
xlabel(ha,'y')
ylabel(ha,'p(y)')
legend(ha,{...
    'prior predictive',...
    'moment-matching',...
    'Laplace',...
    ['split-Laplace: K=',num2str(nmax)]...
    })
ha(2) = subplot(2,2,3,'parent',hf);
set(ha(2),'nextplot','add')
plot(ha(2),muy,'ko')
plot(ha(2),[1 nmog],[Ey0 Ey0],'r')
grid(ha(2),'on')
ylabel(ha(2),'E[y|K]')
xlabel(ha(2),'K')
set(ha(2),'xtick',[1:10])
ha(3) = subplot(2,2,4,'parent',hf);
set(ha(3),'nextplot','add')
plot(ha(3),Vy,'ko')
plot(ha(3),[1 nmog],[Vy0 Vy0],'r')
grid(ha(3),'on')
set(ha(3),'xtick',[1:10])
ylabel(ha(3),'V[y|K]')
xlabel(ha(3),'K')


% fit exponential convergent process to E[y] and V[y]
opt.inG.x = 1:nmax;
opt.inG.up = 0;
opt.priors.muPhi = [0;1;0];
opt.priors.SigmaPhi = 1e4.*eye(3);
g_fname = @g_expConv;
dim = struct('n',0,'n_theta',0,'n_phi',3);
[p1,o1] = VBA_NLStateSpaceModel(muy(1:nmax),[],[],g_fname,dim,opt);
plot(ha(2),[1:nmax],o1.suffStat.gx,'k--')
opt.inG.up = 1;
[p2,o2] = VBA_NLStateSpaceModel(Vy(1:nmax),[],[],g_fname,dim,opt);
plot(ha(3),[1:nmax],o2.suffStat.gx,'k--')

