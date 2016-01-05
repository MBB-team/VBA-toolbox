% this demo exemplifies the use of mixed-effects analysis in VBA

clear all
close all
clc

ns = 8; % #subjects
dim.n_phi = 1;
dim.n = 1;
dim.n_theta = 1;
dim.p = 1; % data dim (within-subject)
dim.n_t = 16; % trials (within-subject)

% simulate MFX
y = cell(ns,1);
nu = ones(3,1); % population mean
alpha = 1; % population precision
sigmas = ones(ns,1); % within-subject residual variance
g_fname = @g_vgo;
f_fname = @f_vgo;
for i=1:ns
    % draw within-subject effects from population distribution
    params(:,i) = nu + randn(3,1)./sqrt(alpha);
    u{i} = randn(1,dim.n_t);
    options{i}.DisplayWin = 0;
    options{i}.verbose = 0;
    options{i}.dim = dim;
    options{i}.binomial = 1;
    [y{i}] = simulateNLSS(dim.n_t,f_fname,g_fname,params(1,i),params(2,i),u{i},Inf,[],options{i},params(3,i));
end


% priors_group.QPhi = 0.*eye(dim.n_phi);
% priors_group.QTheta = 0.*eye(dim.n_theta);
% priors_group.QX0 = 0.*eye(dim.n);
% priors_group.QPhi(2,2) = 0; % ffx
% priors_group.SigmaPhi = eye(dim.n_phi);
% priors_group.SigmaPhi(3,3) = 0; % fix population mean to 0
[p_sub,o_sub,p_group,o_group] = VBA_MFX(y,u,f_fname,g_fname,dim,options);%,priors_group);


