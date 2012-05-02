function [y,post,out,post123,out123] = demo_fool_dcmWhrf(lo,i,j)

% Demo for sDCM for fMRI (with 'missing region')

close all
% clear variables

%-----------------------------------------------------------
%-------------- DCM model specification --------------------

%--- Basic settings
n_t = 0.5e2;                    % number of time samples
TR = 3e0;                       % sampling period (in sec)
microDT = 1e-1;                 % micro-time resolution (in sec)
u       = zeros(2,n_t);         % deterministic inputs
u(1,1:5) = 1;
u(1,30:35) = 1;
u(2,10:15) = 1;
u(2,25:30) = 1;

%--- DCM structure
nreg = 4; % nb regions
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

f_fname = @f_DCMwHRF;
g_fname = @g_HRF3;
[options] = prepare_fullDCM(A,B,C,D,TR,microDT);
dim.n_theta         = options.inF.ind5(end);
dim.n_phi           = 2*nreg;
dim.n               = 5*nreg;


% Simulate DCM with 4 regions
% A matrix
t_A = exp([ -1.5
            -1.5
            -0.5
            -2.5
            lo.theta(j) ]);
% self-inhibition gain
t_Aself = -1;
% B matrices
t_B{1} = [];
t_B{2} = [];
% C matrix
t_C = exp([ +0.1; +0.1 ]);
% D matrices
t_D{1} = [];
t_D{2} = [];
t_D{3} = [];
t_D{4} = [];

%--- simulated evolution parameters: hemodynamic level
t_E0 = zeros(nreg,1);       % HbO2 extraction fraction gain
t_tau0 = zeros(nreg,1);     % mean blood transit time gain
t_kas = zeros(nreg,1);      % vasodilatory signal decay gain
t_kaf = zeros(nreg,1);      % vasodilatory signal feedback regulation
t_alpha = zeros(nreg,1);    % vessel stifness gain

%--- simulated observation parameters
p_epsilon = zeros(nreg,1);  % ratio of intra- and extravascular signal
p_E0 = zeros(nreg,1);       % HbO2 extraction fraction gain

% precision hyperparameters
alpha   = 1e2;
sigma   = lo.sigma(i);

%--- Recollect paramters for simulated data
nu = size(u,1);
theta = zeros(dim.n_theta,1);
theta(options.inF.indA) = t_A;
for i=1:nu
    theta(options.inF.indB{i}) = t_B{i};
end
theta(options.inF.indC) = t_C;
for i=1:nreg
    theta(options.inF.indD{i}) = t_D{i};
end
theta(options.inF.indself) = t_Aself;
theta(options.inF.ind1) = t_E0;
theta(options.inF.ind2) = t_tau0;
theta(options.inF.ind3) = t_kaf;
theta(options.inF.ind4) = t_kas;
theta(options.inF.ind5) = t_alpha;
phi = zeros(dim.n_phi,1);
phi(options.inG.ind1) = p_E0;
phi(options.inG.ind2) = p_epsilon;



% Build priors for model inversion
priors.muX0 = kron(ones(nreg,1),[0;0;0;0;0]);
priors.SigmaX0 = 0e-1*eye(5*nreg);
priors.muTheta = 0*ones(dim.n_theta,1);
priors.muTheta(options.inF.indself) = -1;
priors.SigmaTheta = 1e-3*eye(dim.n_theta);
priors.SigmaTheta(options.inF.ind1,options.inF.ind1) = 0;
priors.muPhi = 0*ones(dim.n_phi,1);
priors.SigmaPhi = 0*1e-3*eye(dim.n_phi);
priors.a_alpha = 1e6;
priors.b_alpha = 1e4;
priors.a_sigma = 1e4;
priors.b_sigma = 1e0;
% Stat noise time-dependent covariance structure.
% NB: ratio of neural versus hemodynamic precision
for t = 1:n_t
    dq = 1e4*ones(dim.n,1);
    dq(options.inF.n5) = 1;
    priors.iQx{t} = diag(dq);
end
options.priors = priors;

% Build options structure for model inversion
options.priors  = priors;
options.DisplayWin = 0;

% Build time series of hidden states and observations
[y,x,x0,eta,e] = simulateNLSS(n_t,f_fname,g_fname,theta,phi,u,alpha,sigma,options);

% display time series of hidden states and observations
% displaySimulations(y,x,eta,e)
% disp('--paused--')
% pause

% Call inversion routine
% u = zeros(size(u));
options.priors = rmfield(options.priors,'iQx');
[post,out] = VBA_NLStateSpaceModel(y,u,f_fname,g_fname,dim,options);

% Display results
% displayResults(post,out,y,x,x0,theta,phi,alpha,sigma)



%----- Now invert sub-DCMs and do model comparison.

% Select data from regions 1,2,3:
ys = y(1:3,:);

% Change input and model specification...
u0 = u;
u = u(1,:);
nreg = 3;
dim.n_phi = 2*nreg;
dim.n = 5*nreg;
B = cell(1,1);
C = [1;0;0];
D = cell(nreg,1);
clear priors;


% Sub-DCM #1:
A = [0 1 0
     1 0 1
     0 1 0];
[options] = prepare_fullDCM(A,B,C,D,TR,microDT);
dim.n_theta         = options.inF.ind5(end);
priors.muX0 = kron(ones(nreg,1),[0;0;0;0;0]);
priors.SigmaX0 = 0e-1*eye(5*nreg);
priors.muTheta = 0*ones(dim.n_theta,1);
priors.muTheta(options.inF.indself) = -1;
priors.SigmaTheta = 1e-3*eye(dim.n_theta);
priors.SigmaTheta(options.inF.ind1,options.inF.ind1) = 0;
priors.muPhi = 0*ones(dim.n_phi,1);
priors.SigmaPhi = 0*1e-3*eye(dim.n_phi);
priors.a_alpha = 1e6;
priors.b_alpha = 1e4;
priors.a_sigma = 1e4;
priors.b_sigma = 1e0;
for t = 1:n_t
    dq = 1e4*ones(dim.n,1);
    dq(options.inF.n5) = 1;
    priors.iQx{t} = diag(dq);
end
options.priors = priors;
options.DisplayWin = 0;
[post123{1},out123{1}] = VBA_NLStateSpaceModel(ys,u,f_fname,g_fname,dim,options);



% Sub-DCM #2:
A = [0 1 0
     1 0 1
     1 1 0];
[options] = prepare_fullDCM(A,B,C,D,TR,microDT);
dim.n_theta         = options.inF.ind5(end);
priors.muX0 = kron(ones(nreg,1),[0;0;0;0;0]);
priors.SigmaX0 = 0e-1*eye(5*nreg);
priors.muTheta = 0*ones(dim.n_theta,1);
priors.muTheta(options.inF.indself) = -1;
priors.SigmaTheta = 1e-3*eye(dim.n_theta);
priors.SigmaTheta(options.inF.ind1,options.inF.ind1) = 0;
priors.muPhi = 0*ones(dim.n_phi,1);
priors.SigmaPhi = 0*1e-3*eye(dim.n_phi);
priors.a_alpha = 1e6;
priors.b_alpha = 1e4;
priors.a_sigma = 1e4;
priors.b_sigma = 1e0;
for t = 1:n_t
    dq = 1e4*ones(dim.n,1);
    dq(options.inF.n5) = 1;
    priors.iQx{t} = diag(dq);
end
options.priors = priors;
options.DisplayWin = 0;
[post123{2},out123{2}] = VBA_NLStateSpaceModel(ys,u,f_fname,g_fname,dim,options);




% Sub-DCM #3:
A = [0 1 1
     1 0 1
     0 1 0];
[options] = prepare_fullDCM(A,B,C,D,TR,microDT);
dim.n_theta         = options.inF.ind5(end);
priors.muX0 = kron(ones(nreg,1),[0;0;0;0;0]);
priors.SigmaX0 = 0e-1*eye(5*nreg);
priors.muTheta = 0*ones(dim.n_theta,1);
priors.muTheta(options.inF.indself) = -1;
priors.SigmaTheta = 1e-3*eye(dim.n_theta);
priors.SigmaTheta(options.inF.ind1,options.inF.ind1) = 0;
priors.muPhi = 0*ones(dim.n_phi,1);
priors.SigmaPhi = 0*1e-3*eye(dim.n_phi);
priors.a_alpha = 1e6;
priors.b_alpha = 1e4;
priors.a_sigma = 1e4;
priors.b_sigma = 1e0;
for t = 1:n_t
    dq = 1e4*ones(dim.n,1);
    dq(options.inF.n5) = 1;
    priors.iQx{t} = diag(dq);
end
options.priors = priors;
options.DisplayWin = 0;
[post123{3},out123{3}] = VBA_NLStateSpaceModel(ys,u,f_fname,g_fname,dim,options);



% Sub-DCM #4:
A = [0 1 1
     1 0 1
     1 1 0];
[options] = prepare_fullDCM(A,B,C,D,TR,microDT);
dim.n_theta         = options.inF.ind5(end);
priors.muX0 = kron(ones(nreg,1),[0;0;0;0;0]);
priors.SigmaX0 = 0e-1*eye(5*nreg);
priors.muTheta = 0*ones(dim.n_theta,1);
priors.muTheta(options.inF.indself) = -1;
priors.SigmaTheta = 1e-3*eye(dim.n_theta);
priors.SigmaTheta(options.inF.ind1,options.inF.ind1) = 0;
priors.muPhi = 0*ones(dim.n_phi,1);
priors.SigmaPhi = 0*1e-3*eye(dim.n_phi);
priors.a_alpha = 1e6;
priors.b_alpha = 1e4;
priors.a_sigma = 1e4;
priors.b_sigma = 1e0;
for t = 1:n_t
    dq = 1e4*ones(dim.n,1);
    dq(options.inF.n5) = 1;
    priors.iQx{t} = diag(dq);
end
options.priors = priors;
options.DisplayWin = 0;
[post123{4},out123{4}] = VBA_NLStateSpaceModel(ys,u,f_fname,g_fname,dim,options);



% Sub-DCM #5:
A = [0 1 0
     1 0 0
     0 1 0];
[options] = prepare_fullDCM(A,B,C,D,TR,microDT);
dim.n_theta         = options.inF.ind5(end);
priors.muX0 = kron(ones(nreg,1),[0;0;0;0;0]);
priors.SigmaX0 = 0e-1*eye(5*nreg);
priors.muTheta = 0*ones(dim.n_theta,1);
priors.muTheta(options.inF.indself) = -1;
priors.SigmaTheta = 1e-3*eye(dim.n_theta);
priors.SigmaTheta(options.inF.ind1,options.inF.ind1) = 0;
priors.muPhi = 0*ones(dim.n_phi,1);
priors.SigmaPhi = 0*1e-3*eye(dim.n_phi);
priors.a_alpha = 1e6;
priors.b_alpha = 1e4;
priors.a_sigma = 1e4;
priors.b_sigma = 1e0;
for t = 1:n_t
    dq = 1e4*ones(dim.n,1);
    dq(options.inF.n5) = 1;
    priors.iQx{t} = diag(dq);
end
options.priors = priors;
options.DisplayWin = 0;
[post123{5},out123{5}] = VBA_NLStateSpaceModel(ys,u,f_fname,g_fname,dim,options);


% Sub-DCM #6:
A = [0 0 0
     1 0 1
     0 1 0];
[options] = prepare_fullDCM(A,B,C,D,TR,microDT);
dim.n_theta         = options.inF.ind5(end);
priors.muX0 = kron(ones(nreg,1),[0;0;0;0;0]);
priors.SigmaX0 = 0e-1*eye(5*nreg);
priors.muTheta = 0*ones(dim.n_theta,1);
priors.muTheta(options.inF.indself) = -1;
priors.SigmaTheta = 1e-3*eye(dim.n_theta);
priors.SigmaTheta(options.inF.ind1,options.inF.ind1) = 0;
priors.muPhi = 0*ones(dim.n_phi,1);
priors.SigmaPhi = 0*1e-3*eye(dim.n_phi);
priors.a_alpha = 1e6;
priors.b_alpha = 1e4;
priors.a_sigma = 1e4;
priors.b_sigma = 1e0;
for t = 1:n_t
    dq = 1e4*ones(dim.n,1);
    dq(options.inF.n5) = 1;
    priors.iQx{t} = diag(dq);
end
options.priors = priors;
options.DisplayWin = 0;
[post123{6},out123{6}] = VBA_NLStateSpaceModel(ys,u,f_fname,g_fname,dim,options);



% Sub-DCM #7:
A = [0 0 0
     1 0 0
     0 1 0];
[options] = prepare_fullDCM(A,B,C,D,TR,microDT);
dim.n_theta         = options.inF.ind5(end);
priors.muX0 = kron(ones(nreg,1),[0;0;0;0;0]);
priors.SigmaX0 = 0e-1*eye(5*nreg);
priors.muTheta = 0*ones(dim.n_theta,1);
priors.muTheta(options.inF.indself) = -1;
priors.SigmaTheta = 1e-3*eye(dim.n_theta);
priors.SigmaTheta(options.inF.ind1,options.inF.ind1) = 0;
priors.muPhi = 0*ones(dim.n_phi,1);
priors.SigmaPhi = 0*1e-3*eye(dim.n_phi);
priors.a_alpha = 1e6;
priors.b_alpha = 1e4;
priors.a_sigma = 1e4;
priors.b_sigma = 1e0;
for t = 1:n_t
    dq = 1e4*ones(dim.n,1);
    dq(options.inF.n5) = 1;
    priors.iQx{t} = diag(dq);
end
options.priors = priors;
options.DisplayWin = 0;
[post123{7},out123{7}] = VBA_NLStateSpaceModel(ys,u,f_fname,g_fname,dim,options);



% Sub-DCM #8:
A = [0 0 0
     0 0 1
     1 0 0];
[options] = prepare_fullDCM(A,B,C,D,TR,microDT);
dim.n_theta         = options.inF.ind5(end);
priors.muX0 = kron(ones(nreg,1),[0;0;0;0;0]);
priors.SigmaX0 = 0e-1*eye(5*nreg);
priors.muTheta = 0*ones(dim.n_theta,1);
priors.muTheta(options.inF.indself) = -1;
priors.SigmaTheta = 1e-3*eye(dim.n_theta);
priors.SigmaTheta(options.inF.ind1,options.inF.ind1) = 0;
priors.muPhi = 0*ones(dim.n_phi,1);
priors.SigmaPhi = 0*1e-3*eye(dim.n_phi);
priors.a_alpha = 1e6;
priors.b_alpha = 1e4;
priors.a_sigma = 1e4;
priors.b_sigma = 1e0;
for t = 1:n_t
    dq = 1e4*ones(dim.n,1);
    dq(options.inF.n5) = 1;
    priors.iQx{t} = diag(dq);
end
options.priors = priors;
options.DisplayWin = 0;
[post123{8},out123{8}] = VBA_NLStateSpaceModel(ys,u,f_fname,g_fname,dim,options);



% Sub-DCM #9:
A = [0 0 1
     0 0 1
     1 0 0];
[options] = prepare_fullDCM(A,B,C,D,TR,microDT);
dim.n_theta         = options.inF.ind5(end);
priors.muX0 = kron(ones(nreg,1),[0;0;0;0;0]);
priors.SigmaX0 = 0e-1*eye(5*nreg);
priors.muTheta = 0*ones(dim.n_theta,1);
priors.muTheta(options.inF.indself) = -1;
priors.SigmaTheta = 1e-3*eye(dim.n_theta);
priors.SigmaTheta(options.inF.ind1,options.inF.ind1) = 0;
priors.muPhi = 0*ones(dim.n_phi,1);
priors.SigmaPhi = 0*1e-3*eye(dim.n_phi);
priors.a_alpha = 1e6;
priors.b_alpha = 1e4;
priors.a_sigma = 1e4;
priors.b_sigma = 1e0;
for t = 1:n_t
    dq = 1e4*ones(dim.n,1);
    dq(options.inF.n5) = 1;
    priors.iQx{t} = diag(dq);
end
options.priors = priors;
options.DisplayWin = 0;
[post123{9},out123{9}] = VBA_NLStateSpaceModel(ys,u,f_fname,g_fname,dim,options);


% Sub-DCM #10:
A = [0 0 0
     0 0 1
     1 1 0];
[options] = prepare_fullDCM(A,B,C,D,TR,microDT);
dim.n_theta         = options.inF.ind5(end);
priors.muX0 = kron(ones(nreg,1),[0;0;0;0;0]);
priors.SigmaX0 = 0e-1*eye(5*nreg);
priors.muTheta = 0*ones(dim.n_theta,1);
priors.muTheta(options.inF.indself) = -1;
priors.SigmaTheta = 1e-3*eye(dim.n_theta);
priors.SigmaTheta(options.inF.ind1,options.inF.ind1) = 0;
priors.muPhi = 0*ones(dim.n_phi,1);
priors.SigmaPhi = 0*1e-3*eye(dim.n_phi);
priors.a_alpha = 1e6;
priors.b_alpha = 1e4;
priors.a_sigma = 1e4;
priors.b_sigma = 1e0;
for t = 1:n_t
    dq = 1e4*ones(dim.n,1);
    dq(options.inF.n5) = 1;
    priors.iQx{t} = diag(dq);
end
options.priors = priors;
options.DisplayWin = 0;
[post123{10},out123{10}] = VBA_NLStateSpaceModel(ys,u,f_fname,g_fname,dim,options);



% Sub-DCM #11:
A = [0 0 1
     0 0 1
     1 1 0];
[options] = prepare_fullDCM(A,B,C,D,TR,microDT);
dim.n_theta         = options.inF.ind5(end);
priors.muX0 = kron(ones(nreg,1),[0;0;0;0;0]);
priors.SigmaX0 = 0e-1*eye(5*nreg);
priors.muTheta = 0*ones(dim.n_theta,1);
priors.muTheta(options.inF.indself) = -1;
priors.SigmaTheta = 1e-3*eye(dim.n_theta);
priors.SigmaTheta(options.inF.ind1,options.inF.ind1) = 0;
priors.muPhi = 0*ones(dim.n_phi,1);
priors.SigmaPhi = 0*1e-3*eye(dim.n_phi);
priors.a_alpha = 1e6;
priors.b_alpha = 1e4;
priors.a_sigma = 1e4;
priors.b_sigma = 1e0;
for t = 1:n_t
    dq = 1e4*ones(dim.n,1);
    dq(options.inF.n5) = 1;
    priors.iQx{t} = diag(dq);
end
options.priors = priors;
options.DisplayWin = 0;
[post123{11},out123{11}] = VBA_NLStateSpaceModel(ys,u,f_fname,g_fname,dim,options);


% % Make predictions
% try
%     options = out.options;
%     [xs,ys,xhat,vx,yhat,vy] = comparePredictions(...
%         n_t,theta,phi,zeros(size(u)),alpha,sigma,options,posterior,dim);
% catch
%     disp('------!!Unable to form predictions!!------')
% end

