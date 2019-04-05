% Demo for DCM for fMRI (distributed responses)
% This demo inverts the DCM for fMRI model, which contains the ballon model
% as a generalized observation function (not affected by stochastic
% innovations).

close all
clear variables

%-----------------------------------------------------------
%-------------- DCM model specification --------------------

%--- Basic settings
TR = 1e0;                      % sampling period (in sec)
n_t = round(1e2/TR);            % number of time samples
dtU = round(2/TR)+1;           % input-on time interval
t0U = round(10/TR)+1;                  
microDT = 1e-1;               % micro-time resolution (in sec)
homogeneous = 1;              % params of g(x) homogeneous accross regions
reduced_f = 1;
alpha   = 2e2/TR;              % state noise precision
sigma   = 1e0;              % measurement noise precision

%--- Input
u       = zeros(2,n_t);
u(1,t0U:t0U+dtU) = 1;
u(1,4*t0U:4*t0U+dtU) = 1;
u(2,4*t0U:4*t0U+5*dtU) = 1;

%--- DCM structure
% invariant effective connectivity
A = [0 1 1
     1 0 1
     0 1 0];
nreg = size(A,1);
% modulatory effects
B{1} = zeros(nreg,nreg);
B{2} = [0 0 0
        1 0 0
        0 0 0];
% input-state coupling
C = [1 0
     0 0
     0 0];
% gating (nonlinear) effects
D{1} = [0 0 0
        0 0 0
        0 1 0];
D{2} = zeros(nreg,nreg);
D{3} = zeros(nreg,nreg);

%--- Build options/dim structures for model inversion
disp('Extracting HRF parameters...')
[thetaHRF,phiHRF]   = get_HRFparams(TR,microDT,homogeneous);
disp('Done.')
thetaHRF = 0.*thetaHRF;
phiHRF = 0.*phiHRF;
f_fname = @f_DCMwHRF;
g_fname = @g_HRF_distributed;
[options] = prepare_fullDCM(A,B,C,D,TR,microDT);
options.inG.n_phi = 4; % number of spatial modes
options.inG.B = VBA_sigmoid(randn(8,options.inG.n_phi)); % spatial modes
options.inG.ind_hrf = 1:2*nreg;
options.inG.n_reg = nreg;
for i=1:nreg
    options.inG.ind3{i} = options.inG.ind_hrf(end)+1+(i-1)*4:...
        options.inG.ind_hrf(end)+i*4;
end
options.MaxIter = 8; % to quicken inversion
dim.n_theta         = options.inF.ind5(end);
dim.n_phi           = (2+options.inG.n_phi).*nreg;
dim.n               = 5*nreg;
dim.p               = 8;
options.dim=dim;

%--- Build priors for model inversion
indHemo = options.inF.indself+1:dim.n_theta;
priors.muX0 = kron(ones(nreg,1),[0;0;0;0;0]);
priors.SigmaX0 = 0e-1*speye(5*nreg);
priors.muTheta = 1e-1*ones(dim.n_theta,1);
priors.muTheta(options.inF.indself) = -0;
priors.SigmaTheta = 1e-0*eye(dim.n_theta);
priors.SigmaTheta(options.inF.indself,options.inF.indself) = 0;
if reduced_f
    % fix some HRF params to their default values
    priors.muTheta(options.inF.ind1) = thetaHRF(1);
    priors.muTheta(options.inF.ind3) = thetaHRF(3);
    priors.muTheta(options.inF.ind5) = thetaHRF(6);
    priors.SigmaTheta(options.inF.ind1,options.inF.ind1) = 0;
    priors.SigmaTheta(options.inF.ind3,options.inF.ind3) = 0;
    priors.SigmaTheta(options.inF.ind5,options.inF.ind5) = 0;
end
priors.muPhi = 1e-1*ones(dim.n_phi,1);
priors.SigmaPhi = 1e-0*eye(dim.n_phi);
priors.SigmaPhi(options.inG.ind_hrf,options.inG.ind_hrf) = 0;

% NB on hyperpriors:
%   - fix state noise precision using high scale param
%   - use non-informative priors on the residual precision, with high
%   expectation.
% This is because of the first iteration of the hidden states posterior
% update, which has to deviate from the its prior predictive density (as
% derived from the deterministic inversion).
% The following iterations will then work with a realistic (expected)
% residual precision, and adapt.

SC = 1e2;
priors.a_alpha = Inf;%SC*alpha;
priors.b_alpha = 0;%SC;
priors.a_sigma = 1e0;
priors.b_sigma = 1e-4;


% State noise time-dependent covariance structure.
% NB: ratio of neural versus hemodynamic precision
for t = 1:n_t
    dq = 1e4*ones(dim.n,1);
    dq(options.inF.n5) = 1;
    priors.iQx{t} = diag(dq);
end
options.priors = priors;
options.backwardLag = 10;
options.gradF = 0;
options.updateHP = 1;
options.init0 = 0;


%-----------------------------------------------------------
%----------- simulated times series specification ----------

%--- simulated evolution parameters: neuronal level
% A matrix
t_A = exp([ -0.5
            -0.5
            -0.5
            -2.5
            -1.5 ]);
% self-inhibition gain
t_Aself = -0;
% B matrices
t_B{1} = [];
t_B{2} = exp([ -0.5 ]);
% C matrix
t_C = exp([ +0.1 ]);
% D matrices
t_D{1} = -exp([ -1.5 ]);
t_D{2} = [];
t_D{3} = [];

%--- simulated evolution parameters: hemodynamic level
t_E0 = thetaHRF(1)*ones(nreg,1);       % HbO2 extraction fraction gain
t_tau0 = thetaHRF(2)*ones(nreg,1);     % mean blood transit time gain
t_kaf = thetaHRF(3)*ones(nreg,1);      % vasodilatory signal feedback regulation
t_kas = thetaHRF(4)*ones(nreg,1);      % vasodilatory signal decay gain
t_alpha = thetaHRF(6)*ones(nreg,1);    % vessel stifness gain

%--- simulated observation parameters
if ~homogeneous
    p_E0 = phiHRF(1)*ones(nreg,1);       % HbO2 extraction fraction gain
    p_epsilon = phiHRF(2)*ones(nreg,1);  % ratio of intra- and extravascular signal
else
    p_E0 = phiHRF(1);
    p_epsilon = phiHRF(2);
end
p_w = repmat([1;2;-1;-2],nreg,1);   % spatial pattern parameters

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
phi(options.inG.ind1) = p_E0;
phi(options.inG.ind2) = p_epsilon;
phi = [ phi(:) ; p_w ];

% Build time series of hidden states and observations
[y,x,x0,eta,e] = VBA_simulate (n_t,f_fname,g_fname,theta,phi,u,alpha,sigma,options);

% display time series of hidden states and observations
displaySimulations(y,x,eta,e);
% disp('--paused--')
% pause

%-----------------------------------------------------------
%------------------- model inversion -----------------------

% Call inversion routine
[posterior,out] = VBA_NLStateSpaceModel(y,u,f_fname,g_fname,dim,options);

% Display results
displayResults(posterior,out,y,x,x0,theta,phi,alpha,sigma);

% Make predictions
try
    options = out.options;
    [xs,ys,xhat,vx,yhat,vy] = VBA_comparePredictions(...
        n_t,theta,phi,u,alpha,sigma,options,posterior,dim);
catch
    disp('------!!Unable to form predictions!!------')
end
