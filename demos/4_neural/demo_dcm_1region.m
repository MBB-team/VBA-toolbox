% Demo for 1-region DCM for fMRI.
% This demo inverts the DCM for fMRI model, which contains the ballon model
% as a generalized observation function (not affected by stochastic
% innovations).

close all
clear variables

%-----------------------------------------------------------
%-------------- DCM model specification --------------------

%--- Basic settings
f_fname = @f_DCMwHRF;
g_fname = @g_HRF3;
TR = 1e0;                      % sampling period (in sec)
n_t = round(6e1/TR);            % number of time samples
dtU = 1;%round(2/TR)+1;           % input-on time interval
t0U = 1;%round(10/TR)+1;                  
microDT = 1e-1;               % micro-time resolution (in sec)
homogeneous = 1;              % params of g(x) homogeneous accross regions
lin = 0;
alpha = 1e2/TR;              % state noise precision
sigma = Inf;%1e2;              % measurement noise precision
stochastic = 0;

%--- Input
u = zeros(1,n_t);%zeros(2,n_t);
times = 2:2+16;
u(1,times) = 1;
% u(2,floor(n_t/2):n_t) = 1;
% % u(1,4*t0U:4*t0U+dtU) = 1;

%--- DCM structure
% invariant effective connectivity
A = [0];
nreg = size(A,1);
% modulatory effects
B{1} = 0;
B{2} = 0;
% input-state coupling
C = 1;%[1,0];
% gating (nonlinear) effects
D{1} = 0;

%--- Build options/dim structures for model inversion
thetaHRF = zeros(6,1); % use SPM prior mean for HRF parameters
phiHRF = zeros(2,1);
[options]           = prepare_fullDCM(A,B,C,D,TR,microDT,homogeneous);
dim.n_theta         = options.inF.ind5(end);
dim.n_phi           = 2;
dim.n               = 5*nreg;

%--- Build priors for model inversion
options.priors = getPriors(nreg,n_t,options,0,stochastic);
options.microU = 0;
options.backwardLag = ceil(16/TR);  % 16 secs effective backward lag
options.inF.linearized = lin;

%-----------------------------------------------------------
%----------- simulated times series specification ----------

%--- simulated evolution parameters: neuronal level
% A matrix
t_A = [];
% self-inhibition gain
t_Aself = -0;
% B matrices
t_B{1} = [];
t_B{2} = [];%-10;
% C matrix
t_C = 1;
% D matrices
t_D{1} = [];

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



%--- Simulate time series of hidden states and observations
[y,x,x0,eta,e] = VBA_simulate (n_t,f_fname,g_fname,theta,phi,u,alpha,sigma,options);

% display time series of hidden states and observations
displaySimulations(y,x,eta,e);
% disp('--paused--')
% pause


%-----------------------------------------------------------
%------------------- model inversion -----------------------
%--- Call inversion routine
[posterior,out] = VBA_NLStateSpaceModel(y,u,f_fname,g_fname,dim,options);

%--- Display results
displayResults(posterior,out,y-e,x,x0,theta,phi,alpha,sigma);

%--- Make predictions
try
    options = out.options;
    [xs,ys,xhat,vx,yhat,vy] = VBA_comparePredictions(...
        n_t,theta,phi,u,alpha,sigma,options,posterior,dim);
catch
    disp('------!!Unable to form predictions!!------')
end


