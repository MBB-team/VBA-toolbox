% Demo: design optimization for DCM of fMRI data.


close all
clear variables



%-----------------------------------------------------------
%-------------- DCM model specification --------------------

%--------- Basic settings ----------------------------------
inG.t0 = 0;
inG.dt = 1e-1;
inG.tf = 1e1;
f_fname = [];
g_fname = @g_sumOfSines;

dim.n = 0;
dim.n_phi = 2;
dim.n_theta = 0;
dim.n_t = 1;
dim.p = inG.tf./inG.dt + 1;

priors.muPhi = zeros(2,1);
priors.SigmaPhi = 1e2*eye(2);
priors0 = priors;
priors0.SigmaPhi(2,2) = 0;

options.inG = inG;
options.dim = dim;
options.priors = priors;
options0 = options;
options0.priors = priors0;

DIM = {dim,dim};
OPT = {options,options0};

gridu = -5:5e-2:5;
for i=1:length(gridu)
    [DJS,out] = designEfficiency(f_fname,g_fname,DIM,OPT,gridu(i),'models');
    bLC(i) = out.b;
end

hf = figure('color',[1 1 1]);
ha = axes('parent',hf);
plot(ha,gridu,bLC,'k.')
xlabel(ha,'phase offset')
ylabel(ha,'Laplace-Chernoff bound')
title(ha,'Comparing m1 and m0:Phi(2)=0: design risk')
box(ha,'off')
set(ha,'ygrid','on')
