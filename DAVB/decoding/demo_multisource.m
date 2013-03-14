function [posterior,out,y] = demo_multisource()
%%% ----------------------------------------------------------------------
%   Initialise the generative model, simulate data and inverse the model
%   
%   IN:     - b ; idyosincratic value of connectivity
%   OUT: classical structure from the DCM toolbox
%%% ----------------------------------------------------------------------

close all

%% -----------------------------------------------------------
% -------------- DCM model specification --------------------

% === Basic settings =======================================
f_fname = @f_DCMwHRFext ; % evolution function (DCM for fMRI + decoding)
g_fname = @g_demo_multisource;  % observation function (multiple sources)

TR = 2e0;                     % sampling period (in sec)
n_t = round(166/TR);          % number of time samples
microDT = 2e-1;               % micro-time resolution (in sec)
homogeneous = 1;              % params of g(x) homogeneous accross regions
reduced_f = 0;                % fix some HRF params
lin = 1;                      % linearized variant of HRF Balloon model
stochastic = 0;               % flag for stochastic DCM inversion
alpha = Inf;                  % state noise precision
sigma = [.5 .5]';%1/2;                  % measurement noise precision

% === Input ================================================
 u = zeros(2,n_t);
% orthogonal modulation
cpt_0=0;
cpt_1=0;
for(i=14:16:n_t)
    if mod(cpt_0,2) <1
        u(1,i:i+8) = .5  ;
    else
        u(1,i:i+8) = .1  ;
    end
    if mod(cpt_1,4) < 2
        u(2,i:i+8) = .5  ;
    else
        u(2,i:i+8) = .1  ;
    end
    cpt_0=cpt_0+1;
    cpt_1=cpt_1+1;
end
% u(:,floor(size(u,2)/2):end)=.5*u(:,floor(size(u,2)/2):end);
nu = size(u,1);

% === DCM structure ========================================
% invariant effective connectivity
A = [0 1;
     1 0];
nreg = size(A,1);
% modulatory effects
B{1} = zeros(nreg,nreg);
B{2} = zeros(nreg,nreg);
% input-state coupling
C = [1 0;
     0 1];
% gating (nonlinear) effects
D{1} = zeros(nreg,nreg);
D{2} = zeros(nreg,nreg);

% === Decoding scheme ========================================
% direct state->response
hA=[1 0;
    0 1];
nrep=size(hA,1);
% modulatory state->response inputs
hB = {zeros(nrep,nreg),zeros(nrep,nreg)}; 
% direct input->response
hC = [];
% gating (nonlinear)
hD = {zeros(nrep,nreg),zeros(nrep,nreg)};

% === Build options and dim structures =====================
options = prepare_fullDCM(A,B,C,D,TR,microDT,homogeneous,hA,hB,hC,hD);
%- sourcing
sources(1) = struct('out',1,'type',0); % BOLD signal (gaussian, dim=2)
sources(2) = struct('out',2,'type',0); % BOLD signal (gaussian, dim=2)
sources(3) = struct('out',3,'type',1); % first binary response
sources(4) = struct('out',4,'type',1); % second binary response
options.sources=sources;
options.priors = getPriors(nreg,n_t,options,reduced_f,stochastic);
options.microU = 0;
options.backwardLag = 8;
options.GnFigs = 0;
options.inF.linearized = lin;
options.DisplayWin=1;
options.verbose=1;

dim.n_theta = options.inF.ind5(end);
if options.extended
    dim.n_phi = options.inG.indr;
else
    dim.n_phi = options.inG.ind2(end);
end
dim.n = 5*nreg+nrep;



%% -----------------------------------------------------------
%----------- simulated times series specification ----------

%--- simulated evolution parameters: neuronal level
theta = zeros(dim.n_theta,1);
phi = zeros(dim.n_phi,1);

%- DCM
t_A=[0 0];
t_Aself = -.1;
t_B{1} = [];
t_B{2} = [] ; 
t_C = [1 1];
t_D{1} = [];
t_D{2} = [];
%- Decoding
t_hA = 2*[1 1]; 
t_hself = log(1/options.inF.deltat) ; % dirac kernel, fixed in priors
t_hB{1} = [];
t_hB{2} = [];
t_hC = [];
t_hD{1} = [];
t_hD{2} = [];

% store evolution parameters
theta(options.inF.indA) = t_A;
theta(options.inF.indself) = t_Aself;
for i=1:nu
    theta(options.inF.indB{i}) = t_B{i};
end
theta(options.inF.indC) = t_C;
for i=1:nreg
    theta(options.inF.indD{i}) = t_D{i};
end
theta(options.inF.indhA) = t_hA;
theta(options.inF.indhself) = t_hself;
for i=1:nu
    theta(options.inF.indhB{i}) = t_hB{i};
end
theta(options.inF.indhC) = t_hC;
for i=1:nreg
    theta(options.inF.indhD{i}) = t_hD{i};
end

%--- simulated observation parameters
phi(options.inG.indr) = 0; % log inverse temperature

%--- Simulate time series of hidden states and observations
disp('*** Simulation');
[y,x,x0,eta,e] = simulateNLSS(n_t,f_fname,g_fname,theta,phi,u,alpha,sigma,options);
% display time series of hidden states and observations
% displaySimulations(y,x,eta,e)
% disp('--paused--')
% pause

%% -----------------------------------------------------------
%  ------------------- model inversion -----------------------

%--- Call inversion routine
disp('*** Hypothesis inversion');

%- data exclusion
options.isYout=zeros(size(y));
options.isYout(1:2,:)=0; % set to 1 to try input-behavior mapping without fMRI !

%- model inversion
[posterior,out] = VBA_NLStateSpaceModel(y,u,f_fname,g_fname,dim,options);

%- 
displayResults(posterior,out,y-e,x,x0,theta,phi,alpha,sigma)


end