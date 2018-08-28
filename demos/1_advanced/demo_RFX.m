% demo for mixed-effects analysis (2-levels hierachical model)
% This script simulates and analyses data under two different hierarchical
% models, namely (i) H0: group-mean=0, and (ii) H1: group-mean~=0.
% The effect of H1 and H0 onto both the group-mean parameters estimates
% (given H1) and the log-Bayes factors (comparisons of H1 vs H0) is
% assessed using Monte-Carlo simulations.
% The trick for performing a mixed-effects analysis using the structure of
% the generative model of VBA is to use hidden states, where "initial
% conditions" (x_0) serve as the group mean, x_1 are subject-dependant
% effects, and the precision (alpha) of the transition density
% p(x_1|x_0,alpha) captures the between-subject variance. This trick
% however, can be used here because we deal with a static model (a
% within-subject GLM). 

clear variables
close all

% Choose basic settings for simulations
ns = 32; % number of subjects
n = 256; % number of observations per subject
f_fname = @g_RFX;
g_fname = @g_RFX;
alpha = 1e1;
sigma = 1e0;
theta = [];
phi = [];
x0 = zeros(ns,1);

X2 = [ones(ns,1),zeros(ns,ns-1)]; % 2d-level design matrix
u = [];

% Build options structure for temporal integration of SDE
inF.ns = ns;
inF.X = X2;
options.inF = inF;

% Build priors for model inversion
priors.muX0 = zeros(ns,1);
priors.SigmaX0 = 0*eye(ns);
priors.a_alpha = 1e0;
priors.b_alpha = 1e0;
priors.a_sigma = 1e0;
priors.b_sigma = 1e0;

% Build options and dim stuctures for model inversion
options.priors = priors;
options.DisplayWin = 0;
dim.n_theta = 0;
dim.n_phi = 0;
dim.n = ns;


Nmcmc = 32;
p = cell(2,2,Nmcmc);
o = cell(2,2,Nmcmc);
F = zeros(2,2,Nmcmc);
X0 = zeros(2,Nmcmc);
for ii=1:Nmcmc
    % simulate data with and without 2nd-level effect
    X1 = randn(n,ns); % 1st-level design matrix
    options.inG.X = X1;
    for i=1:2
        x0(1) = 2-i;
        [y,x,x0,eta,e] = VBA_simulate (1,f_fname,g_fname,theta,phi,u,alpha,sigma,options,x0);
        % Invert model with and without 2nd-level effect
        for j=1:2
            options.priors.SigmaX0(1,1) = 2-j; % group mean effect
            [p{i,j,ii},o{i,j,ii}] = VBA_NLStateSpaceModel(y,u,f_fname,g_fname,dim,options);
            F(i,j,ii) = o{i,j,ii}.F;
            if j==1
                % extract estimaetd group-mean
                X0(i,ii) = p{i,j,ii}.muX0(1);
            end
        end
    end
end


hf = figure('color',[1 1 1]);
pos = get(hf,'position');
set(hf,'position',pos.*[1 1 1.5 1]);
models = {'H1','H0'};

ha = subplot(1,2,1,'parent',hf);
mx0 = mean(X0,2);
vx0 = var(X0,[],2)./Nmcmc;
plotUncertainTimeSeries(mx0,vx0,[],ha);
set(ha,'xlim',[0,3],'xtick',[1,2],'xticklabels',models)
xlabel(ha,'type of simulated data')
ylabel(ha,['E[X0|y,H1]'])
title(ha,'estimated group effect (under H1)')
box(ha,'off')


ha = subplot(1,2,2,'parent',hf);
dF = F(:,1,:) - F(:,2,:);
mdF = mean(dF,3);
vdF = var(dF,[],3)./Nmcmc;
plotUncertainTimeSeries(mdF,vdF,[],ha);
set(ha,'xlim',[0,3],'xtick',[1,2],'xticklabels',models)
xlabel(ha,'type of simulated data')
ylabel(ha,['log p(y|',models{1},') - log p(y|',models{2},')'])
title(ha,'evidence for a group effect')
box(ha,'off')

VBA_getSubplots ();


