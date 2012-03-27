% Demo for Gaussian convolution observation model.
%------------------------------------------------------------
% Copyright (C) 2012 Jean Daunizeau / License GNU GPL v2
%------------------------------------------------------------

close all
clear variables

%---- simulate noisy gbf ----%

% Choose basic settings for simulations
nmax = 10;
p = 200;
sigma = 1;            % precision
phi = [1];         % observation parameters

P0 = -0e-1;
V0 = 1e0;

g_fname = @g_x2; % observation function
inG.x = ones(p,1);  % grid on which the gbf is evaluated


% Build simulated observations
[gx] = feval(g_fname,[],phi,[],inG);
y = gx + sqrt(sigma.^-1)*randn(size(gx));

% display time series of hidden states and observations
hf = figure('color',[1 1 1]);
dp = 1e-3;
gridP = [-4:dp:4];
np = length(gridP);
py = zeros(np,1);
pp0 = py;
pp = py;
for i=1:np
    ey = y-feval(g_fname,[],gridP(i),[],inG);
    py(i) = -0.5.*sigma.*ey'*ey;
    ephi = gridP(i) - P0;
    pp0(i) = -0.5*ephi.^2/V0;
    pp(i) = py(i) + pp0(i);
end
py = exp(py).*(2*pi/sigma)^(-p/2);
pp0 = exp(pp0).*(2*pi*V0)^(-1/2);
pp = exp(pp).*((2*pi/sigma)^(-p/2)).*((2*pi*V0)^(-1/2));
sumpp = sum(pp);
pp = pp./sumpp;
Epp = sum(gridP(:).*pp(:));
Vpp = sum((gridP(:)-Epp).^2.*pp(:));
ppmm = exp(-0.5.*(gridP-Epp).^2./Vpp).*(2*pi*Vpp)^(-1/2);
ha = subplot(2,1,1,'parent',hf);
set(ha,'nextplot','add')
plot(ha,gridP,py./sum(py),'g')
plot(ha,gridP,pp0./sum(pp0),'b')
plot(ha,gridP,pp./sum(pp),'r')
plot(ha,gridP,ppmm./sum(ppmm),'r--')
xlabel(ha,'parameter Phi')
ylabel(ha,'q(Phi)')
% title(ha,'information')
legend(ha,...
    {...
    'likelihood',...
    'prior',...
    'posterior',...
    'moment-matching'})
axis(ha,'tight')
grid(ha,'on')
logEv = log(sumpp) + log(dp);


%---- Invert gbf on simulated data ----%

% Build priors structure
priors.muPhi = P0;
priors.SigmaPhi = V0;
priors.a_sigma = sigma;
priors.b_sigma = 1;
% Build options structure
% options.checkGrads = 1;
options.priors = priors; 
options.inG = inG; 
options.GnFigs = 0;
options.verbose = 1;
options.updateHP = 0;
dim.n_phi = 1;
dim.n_theta = 0;
dim.n=0;
% Call inversion routine
post = cell(nmax,1);
out = cell(nmax,1);
ep = zeros(nmax,1);
vp = zeros(nmax,1);
F = zeros(nmax,1);
for nmog = 1:nmax
    nmog
    options.nmog = nmog;
    [post{nmog},out{nmog}] = ...
        VBA_NLStateSpaceModel(y,[],[],g_fname,dim,options);
    ep(nmog) = post{nmog}.muPhi;
    vp(nmog) = post{nmog}.SigmaPhi;
    F(nmog) = out{nmog}.F;
end

plap = -0.5.*(gridP(:)-post{1}.muPhi).^2./post{1}.SigmaPhi;
plap = exp(plap-max(plap));
plap = plap./sum(plap);

pmog = -0.5.*(gridP(:)-post{nmax}.muPhi).^2./post{nmax}.SigmaPhi;
pmog = exp(pmog-max(pmog));
pmog = pmog./sum(pmog);

plot(ha,gridP,plap,'k:')

plot(ha,gridP,pmog,'k--')
legend(ha,...
    {...
    'likelihood',...
    'prior',...
    'posterior',...
    'moment-matching',...
    'VBL',...
    ['VBsL :K=',num2str(nmax)]})
axis(ha,'tight')
grid(ha,'on')


ha(2) = subplot(2,2,3,'parent',hf);
set(ha(2),'nextplot','add')
plot(ha(2),ep,'ko')
plot(ha(2),[1 nmog],[Epp Epp],'r')
grid(ha(2),'on')
ylabel(ha(2),'E[Phi|y,K]')
xlabel(ha(2),'K')
set(ha(2),'xtick',[1:10])
ha(3) = subplot(2,2,4,'parent',hf);
set(ha(3),'nextplot','add')
plot(ha(3),vp,'ko')
plot(ha(3),[1 nmog],[Vpp Vpp],'r')
grid(ha(3),'on')
set(ha(3),'xtick',[1:10])
ylabel(ha(3),'V[Phi|y,K]')
xlabel(ha(3),'K')


%---- Display results ----%
displayResults(post{nmax},out{nmax},y,[],[],[],phi,[],sigma)

