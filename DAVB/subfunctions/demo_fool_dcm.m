function [y,post,out,post123,out123] = demo_fool_dcm(lo,i,j)

% Demo for sDCM for fMRI (with 'missing region')
% This demo inverts the DCM for fMRI model, without the ballon model, which
% is replaced by a nonlinear sigmoid observation function.

% warning on
% warning('This requires a .mat file containing simulations parameters!')
% return

close all
% clear variables

% Simulate DCM with 4 regions

n_t = 5e2;
f_fname = @f_dcm4fmri;
g_fname = @g_Id;
nreg = 4; % nb regions
nu = 2; % nb inputs
decim = 4;
deltat = 1e0./decim;
u = zeros(2,n_t);
u(1,100:150) = 1;
u(1,300:350) = 1;
u(2,50:100) = 1;
u(2,200:250) = 1;
alpha   = 1e1;%lo.alpha(i);
sigma   = lo.sigma(i);
theta   = [ exp([ -0.5
                  -2
                  -0.5
                  -1
                  lo.theta(j)%-0.5
                  .1
                  .1])
            -2.5        ];
phi     = [];

% DCM specification
A = [0 1 0 1
     1 0 1 0
     0 1 0 0
     0 0 0 0];
B = cell(2,1);
C = [1 0
     0 0
     0 0
     0 1];
D = cell(nreg,1);


% Get model dimensions
dim.n_theta = length(theta);
dim.n_phi = 0;
dim.n = nreg;

% Prepare precalculated matrices for DCM inversion
[inF] = prepare_dcm(A,B,C,D);
inF.deltat = deltat;
inG.k = 1;


% Build priors for model inversion
priors.muX0 = 0*ones(nreg,1);
priors.SigmaX0 = 1e-3*eye(nreg);
priors.muTheta = 0*ones(size(theta));
priors.SigmaTheta = 1e-3*eye(length(theta));
priors.a_alpha = 1e3;
priors.b_alpha = 1e0;
priors.a_sigma = 1e3;
priors.b_sigma = 1e0;
for t = 1:n_t
    dq = 1e4*ones(dim.n,1);
    dq(4) = 1;
    priors.iQx{t} = diag(dq);
end

% Build options structure for model inversion
options.priors  = priors;
options.inF     = inF;
options.inG     = inG;
options.decim   = decim;
options.DisplayWin = 0;

% Build time series of hidden states and observations
[y,x,x0,eta,e] = simulateNLSS(n_t,f_fname,g_fname,theta,phi,u,alpha,sigma,options);

% display time series of hidden states and observations
displaySimulations(y,x,eta,e)
% disp('--paused--')
% pause

% Call inversion routine
% u = zeros(size(u));
[post,out] = VBA_NLStateSpaceModel(y,u,f_fname,g_fname,dim,options);

% Display results
displayResults(post,out,y,x,x0,theta,phi,alpha,sigma)



%----- Now invert sub-DCMs and do model comparison.

% Select data from regions 1,2,3:
[gg,dg] = g_Id(zeros(nreg,1),[],[],inG);
ind123 = find(sum(dg(1:3,:),1));
ys = y(ind123,:);

% Change input and priors
u0 = u;
u = u(1,:);
options.priors = rmfield(options.priors,'iQx');

% Change options structure
options.priors.muX0         = 0*ones(nreg-1,1);
options.priors.SigmaX0      = 1e-3*eye(nreg-1);
options.priors.a_alpha = 1e3;
options.priors.b_alpha = 1e0;
B = cell(1,1);
C = [1;0;0];
D = cell(nreg,1);
dim.n = 3;

% Sub-DCM #1:
A = [0 1 0
     1 0 1
     0 1 0];
[inF]                       = prepare_dcm(A,B,C,D);
inF.deltat                  = deltat;
options.inF                 = inF;
options.priors.muTheta      =0*ones(inF.indself,1);
options.priors.SigmaTheta   = 1e-3*eye(inF.indself);
dim.n_theta                 = inF.indself;
[post123{1},out123{1}] = VBA_NLStateSpaceModel(ys,u,f_fname,g_fname,dim,options);



% Sub-DCM #2:
A = [0 1 0
     1 0 1
     1 1 0];
[inF]                       = prepare_dcm(A,B,C,D);
inF.deltat                  = deltat;
options.inF                 = inF;
options.priors.muTheta      = 0*ones(inF.indself,1);
options.priors.SigmaTheta   = 1e-3*eye(inF.indself);
dim.n_theta                 = inF.indself;
[post123{2},out123{2}]      = ...
    VBA_NLStateSpaceModel(ys,u,f_fname,g_fname,dim,options);



% Sub-DCM #3:
A = [0 1 1
     1 0 1
     0 1 0];
[inF]                       = prepare_dcm(A,B,C,D);
inF.deltat                  = deltat;
options.inF                 = inF;
options.priors.muTheta      = 0*ones(inF.indself,1);
options.priors.SigmaTheta   = 1e-3*eye(inF.indself);
dim.n_theta                 = inF.indself;
[post123{3},out123{3}]      = ...
    VBA_NLStateSpaceModel(ys,u,f_fname,g_fname,dim,options);


% Sub-DCM #4:
A = [0 1 1
     1 0 1
     1 1 0];
[inF]                       = prepare_dcm(A,B,C,D);
inF.deltat                  = deltat;
options.inF                 = inF;
options.priors.muTheta      = 0*ones(inF.indself,1);
options.priors.SigmaTheta   = 1e-3*eye(inF.indself);
dim.n_theta                 = inF.indself;
[post123{4},out123{4}]      = ...
    VBA_NLStateSpaceModel(ys,u,f_fname,g_fname,dim,options);


% Sub-DCM #5:
A = [0 1 0
     1 0 0
     0 1 0];
[inF]                       = prepare_dcm(A,B,C,D);
inF.deltat                  = deltat;
options.inF                 = inF;
options.priors.muTheta      = 0*ones(inF.indself,1);
options.priors.SigmaTheta   = 1e-3*eye(inF.indself);
dim.n_theta                 = inF.indself;
[post123{5},out123{5}]      = ...
    VBA_NLStateSpaceModel(ys,u,f_fname,g_fname,dim,options);

% Sub-DCM #6:
A = [0 0 0
     1 0 1
     0 1 0];
[inF]                       = prepare_dcm(A,B,C,D);
inF.deltat                  = deltat;
options.inF                 = inF;
options.priors.muTheta      = 0*ones(inF.indself,1);
options.priors.SigmaTheta   = 1e-3*eye(inF.indself);
dim.n_theta                 = inF.indself;
[post123{6},out123{6}]      = ...
    VBA_NLStateSpaceModel(ys,u,f_fname,g_fname,dim,options);


% Sub-DCM #7:
A = [0 0 0
     1 0 0
     0 1 0];
[inF]                       = prepare_dcm(A,B,C,D);
inF.deltat                  = deltat;
options.inF                 = inF;
options.priors.muTheta      = 0*ones(inF.indself,1);
options.priors.SigmaTheta   = 1e-3*eye(inF.indself);
dim.n_theta                 = inF.indself;
[post123{7},out123{7}]      = ...
    VBA_NLStateSpaceModel(ys,u,f_fname,g_fname,dim,options);


% Sub-DCM #8:
A = [0 0 0
     0 0 1
     1 0 0];
[inF]                       = prepare_dcm(A,B,C,D);
inF.deltat                  = deltat;
options.inF                 = inF;
options.priors.muTheta      = 0*ones(inF.indself,1);
options.priors.SigmaTheta   = 1e-3*eye(inF.indself);
dim.n_theta                 = inF.indself;
[post123{8},out123{8}]      = ...
    VBA_NLStateSpaceModel(ys,u,f_fname,g_fname,dim,options);


% Sub-DCM #9:
A = [0 0 1
     0 0 1
     1 0 0];
[inF]                       = prepare_dcm(A,B,C,D);
inF.deltat                  = deltat;
options.inF                 = inF;
options.priors.muTheta      = 0*ones(inF.indself,1);
options.priors.SigmaTheta   = 1e-3*eye(inF.indself);
dim.n_theta                 = inF.indself;
[post123{9},out123{9}]      = ...
    VBA_NLStateSpaceModel(ys,u,f_fname,g_fname,dim,options);

% Sub-DCM #10:
A = [0 0 0
     0 0 1
     1 1 0];
[inF]                       = prepare_dcm(A,B,C,D);
inF.deltat                  = deltat;
options.inF                 = inF;
options.priors.muTheta      = 0*ones(inF.indself,1);
options.priors.SigmaTheta   = 1e-3*eye(inF.indself);
dim.n_theta                 = inF.indself;
[post123{10},out123{10}]      = ...
    VBA_NLStateSpaceModel(ys,u,f_fname,g_fname,dim,options);


% Sub-DCM #11:
A = [0 0 1
     0 0 1
     1 1 0];
[inF]                       = prepare_dcm(A,B,C,D);
inF.deltat                  = deltat;
options.inF                 = inF;
options.priors.muTheta      = 0*ones(inF.indself,1);
options.priors.SigmaTheta   = 1e-3*eye(inF.indself);
dim.n_theta                 = inF.indself;
[post123{11},out123{11}]      = ...
    VBA_NLStateSpaceModel(ys,u,f_fname,g_fname,dim,options);

% % Make predictions
% try
%     options = out.options;
%     [xs,ys,xhat,vx,yhat,vy] = comparePredictions(...
%         n_t,theta,phi,zeros(size(u)),alpha,sigma,options,posterior,dim);
% catch
%     disp('------!!Unable to form predictions!!------')
% end

