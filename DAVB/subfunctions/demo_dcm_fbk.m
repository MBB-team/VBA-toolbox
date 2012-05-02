% Demo for DCM for fMRI.
% This demo operates a 2X2 factorial analysis of the influence of feedback
% conections onto a two-nodes network:
% - generates data with and without feedback
% - inverts model with and without feedback
% - perform model comparison and inspect data covariance matrix

close all
clear variables


%--- Basic settings
f_fname = @f_DCMwHRF;
g_fname = @g_HRF3;
TR = 1e0;                      % sampling period (in sec)
n_t = round(3e2/TR);            % number of time samples
microDT = 1e-1;               % micro-time resolution (in sec)
homogeneous = 1;              % params of g(x) homogeneous accross regions
reduced_f = 1;
lin = 0;
alpha = Inf;%2e2/TR;              % state noise precision
sigma = 1e2;              % measurement noise precision
thetaHRF = zeros(6,1);
phiHRF = zeros(2,1);
nreg = 2;

%--- Input
u = sample_u(1,16,0,1,ceil(n_t./TR),16,4);
nu = size(u,1);



%-----------------------
%- model 1: feedback-

A = [0 0
    1 0];
B{1} = zeros(nreg,nreg);
C = [1
    0];
D{1} = zeros(nreg,nreg);
D{2} = zeros(nreg,nreg);
[options{1},dim{1}] = getOptions4dcm(...
    A,...
    B,...
    C,...
    D,...
    TR,...
    microDT,...
    n_t,...
    homogeneous,...
    reduced_f,...
    lin);
options{1}.priors.muTheta = 0.*options{1}.priors.muTheta;
options{1}.priors.SigmaTheta(options{1}.inF.indself,options{1}.inF.indself) = 1e-1;
% options{1}.priors.a_sigma = 1e0;
% options{1}.priors.b_sigma = 1e-4;

t_A = 0.5;
t_Aself = -0;
t_B{1} = [];
t_C = 1;
t_D{1} = [];
t_D{2} = [];
t_E0 = thetaHRF(1)*ones(nreg,1);       % HbO2 extraction fraction gain
t_tau0 = thetaHRF(2)*ones(nreg,1);     % mean blood transit time gain
t_kaf = thetaHRF(3)*ones(nreg,1);      % vasodilatory signal feedback regulation
t_kas = thetaHRF(4)*ones(nreg,1);      % vasodilatory signal decay gain
t_alpha = thetaHRF(6)*ones(nreg,1);    % vessel stifness gain
if ~homogeneous
    p_E0 = phiHRF(1)*ones(nreg,1);       % HbO2 extraction fraction gain
    p_epsilon = phiHRF(2)*ones(nreg,1);  % ratio of intra- and extravascular signal
else
    p_E0 = phiHRF(1);
    p_epsilon = phiHRF(2);
end
theta = zeros(dim{1}.n_theta,1);
theta(options{1}.inF.indA) = t_A;
for i=1:nu
    theta(options{1}.inF.indB{i}) = t_B{i};
end
theta(options{1}.inF.indC) = t_C;
for i=1:nreg
    theta(options{1}.inF.indD{i}) = t_D{i};
end
theta(options{1}.inF.indself) = t_Aself;
theta(options{1}.inF.ind1) = t_E0;
theta(options{1}.inF.ind2) = t_tau0;
theta(options{1}.inF.ind3) = t_kaf;
theta(options{1}.inF.ind4) = t_kas;
theta(options{1}.inF.ind5) = t_alpha;
phi = zeros(dim{1}.n_phi,1);
phi(options{1}.inG.ind1) = p_E0;
phi(options{1}.inG.ind2) = p_epsilon;


[y{1},x{1},x0{1},eta{1},e{1}] = simulateNLSS(...
    n_t,f_fname,g_fname,theta,phi,u,alpha,sigma,options{1});

displaySimulations(y{1},x{1},eta{1},e{1})
Theta{1} = theta;
Phi{1} = phi;
% pause


%-----------------------
%- model 2: feedback+

A = [0 1
    1 0];
B{1} = zeros(nreg,nreg);
C = [1
    0];
D{1} = zeros(nreg,nreg);
D{2} = zeros(nreg,nreg);
[options{2},dim{2}] = getOptions4dcm(...
    A,...
    B,...
    C,...
    D,...
    TR,...
    microDT,...
    n_t,...
    homogeneous,...
    reduced_f,...
    lin);
options{2}.priors.muTheta = 0.*options{2}.priors.muTheta;
options{2}.priors.SigmaTheta(options{2}.inF.indself,options{2}.inF.indself) = 1e-1;
% options{1}.priors.a_sigma = 1e0;
% options{1}.priors.b_sigma = 1e-4;

t_A = [0.5;0.5];
t_Aself = -0;
t_B{1} = [];
t_C = 1;
t_D{1} = [];
t_D{2} = [];
t_E0 = thetaHRF(1)*ones(nreg,1);       % HbO2 extraction fraction gain
t_tau0 = thetaHRF(2)*ones(nreg,1);     % mean blood transit time gain
t_kaf = thetaHRF(3)*ones(nreg,1);      % vasodilatory signal feedback regulation
t_kas = thetaHRF(4)*ones(nreg,1);      % vasodilatory signal decay gain
t_alpha = thetaHRF(6)*ones(nreg,1);    % vessel stifness gain
if ~homogeneous
    p_E0 = phiHRF(1)*ones(nreg,1);       % HbO2 extraction fraction gain
    p_epsilon = phiHRF(2)*ones(nreg,1);  % ratio of intra- and extravascular signal
else
    p_E0 = phiHRF(1);
    p_epsilon = phiHRF(2);
end
theta = zeros(dim{2}.n_theta,1);
theta(options{2}.inF.indA) = t_A;
for i=1:nu
    theta(options{2}.inF.indB{i}) = t_B{i};
end
theta(options{2}.inF.indC) = t_C;
for i=1:nreg
    theta(options{2}.inF.indD{i}) = t_D{i};
end
theta(options{2}.inF.indself) = t_Aself;
theta(options{2}.inF.ind1) = t_E0;
theta(options{2}.inF.ind2) = t_tau0;
theta(options{2}.inF.ind3) = t_kaf;
theta(options{2}.inF.ind4) = t_kas;
theta(options{2}.inF.ind5) = t_alpha;
phi = zeros(dim{2}.n_phi,1);
phi(options{2}.inG.ind1) = p_E0;
phi(options{2}.inG.ind2) = p_epsilon;


[y{2},x{2},x0{2},eta{2},e{2}] = simulateNLSS(...
    n_t,f_fname,g_fname,theta,phi,u,alpha,sigma,options{2});

displaySimulations(y{2},x{2},eta{2},e{2})
Theta{2} = theta;
Phi{2} = phi;
% pause


%-----------------------------------------------------------
%------------------- model inversion -----------------------

for i=1:2
    for j=1:2
        
        [posterior{i,j},out{i,j}] = VBA_NLStateSpaceModel(...
            y{i},...
            u,...
            f_fname,...
            g_fname,...
            dim{j},...
            options{j});
        
        options0 = options{j};
        options0.priors = posterior{i,j};
        
        [muy{i,j},Vy{i,j}] = getLaplace(...
            u,...
            f_fname,...
            g_fname,...
            dim{j},...
            options0,...
            0);
        
        F(i,j) = out{i,j}.F;
        
%         %--- Display results
%         displayResults(...
%             posterior{i,j},...
%             out{i,j},...
%             y{i},...
%             x{i},...
%             x0{i},...
%             Theta{i},...
%             Phi{i},...
%             alpha,...
%             sigma)
        
    end
end

% plot 2X2 model comparison using log Bayes factors
hf = figure('color',[1 1 1]);
ha = axes('parent',hf);
hb = bar(ha,diff(F,[],2));
set(hb,'barwidth',0.6,'facecolor',0.5*[1 1 1])
set(ha,'xticklabel',{'fbk-','fbk+'})
ylabel('log p(y|fbk+) - log p(y|fbk-)')
xlabel('datasets')
set(ha,'box','off','ygrid','on');




su = ~~sum(u,1);
ind = find(diff(su)==-1);
ind2 = find(diff(su)==1);
begin = ind2(ind2>ind(1));
between = begin(1)-1;

begin = max([1,2*(ind2(1)-4)]);
stop = 2*(ind2(2)-1);

C2 = cov2corr(Vy{2,2}(begin:stop,begin:stop));
C1 = cov2corr(Vy{1,1}(begin:stop,begin:stop));

hf = figure ('color',[1 1 1]);

ha(1) = subplot(2,2,1,'parent',hf,'nextplot','add');
dcm_plot_reorderedCov(C1,2,u(ind2(1)-2:end),ha(1));
title(ha(1),'corr(y|fbk-)')

ha(2) = subplot(2,2,2,'parent',hf,'nextplot','add');
dcm_plot_reorderedCov(C2,2,u(ind2(1)-2:end),ha(2));
title(ha(2),'corr(y|fbk+)')

ha(3) = subplot(2,2,3,'parent',hf,'nextplot','add');
dcm_plot_reorderedCov(C2-C1,2,u(ind2(1)-2:end),ha(3));
title(ha(3),'corr(y|fbk+) - corr(y|fbk-)')
getSubplots










