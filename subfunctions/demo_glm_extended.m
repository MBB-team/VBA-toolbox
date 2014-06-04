function [posterior,out] = demo_glm_extended()
%%% ----------------------------------------------------------------------
%   Initialise the generative model, simulate data and inverse the model
%   
%%% ----------------------------------------------------------------------

close all

%-----------------------------------------------------------
%-------------- DCM model specification --------------------

% === Basic settings =======================================
f_fname = @f_DCMwHRFext ;  %
g_fname = @g_demo_extended;   % 

                              % 
TR = 1e0;                     % sampling period (in sec)
n_t = round(64/TR);          % number of time samples
microDT = 1e-1;               % micro-time resolution (in sec)
homogeneous = 1;              % params of g(x) homogeneous accross regions
reduced_f = 0;                % fix some HRF params
lin = 1;                      % linearized variant of HRF Balloon model
stochastic = 0;               % flag for stochastic DCM inversion
alpha = Inf;                  % state noise precision
sigma = Inf;                  % measurement noise precision

% === Input ================================================


u = zeros(1,30);
cpt_0=0;
cpt_1=0;

u(15)=1;
u=repmat(u,1,2);

nu = size(u,1);
n_t = size(u,2);


% === DCM structure ========================================
% invariant effective connectivity

A = [0];
  
nreg = size(A,1);
% modulatory effects
B{1} = zeros(nreg,nreg);
% input-state coupling
C = [1];
% gating (nonlinear) effects
D{1} = zeros(nreg,nreg);

% === Decoding scheme ========================================
hA = [1];


nrep=size(hA,1);
hB = {zeros(nrep,nreg)}; 
hC = zeros(nrep,nu);
hD = {zeros(nrep,nreg)};

% === Build options and dim structures =====================
options = prepare_fullDCM(A,B,C,D,TR,microDT,homogeneous,hA,hB,hC,hD);

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

sources(1) = struct('out',1,'type',0); % BOLD signal (gaussian, dim=4)
sources(2) = struct('out',2,'type',1); % first binary response
options.sources=sources;


%% -----------------------------------------------------------
%----------- simulated times series specification ----------

%--- simulated evolution parameters: neuronal level
theta = zeros(dim.n_theta,1);
phi = zeros(dim.n_phi,1);

%- DCM
t_Aself = .5;
t_A = [0];
t_B{1} = []; 
t_C = [1];
t_D{1} = [];

theta(options.inF.indA) = t_A;
theta(options.inF.indself) = t_Aself;
for i=1:nu
    theta(options.inF.indB{i}) = t_B{i};
end
theta(options.inF.indC) = t_C;
for i=1:nreg
    theta(options.inF.indD{i}) = t_D{i};
end
%- decoding 

t_hA = [1];
t_hB{1} = [];
t_hC = [];
t_hD{1} = [];
t_hself = log(1/options.inF.deltat) ;

theta(options.inF.indhA) = t_hA;
theta(options.inF.indhself) = t_hself;
for i=1:nu
    theta(options.inF.indhB{i}) = t_hB{i};
end
theta(options.inF.indhC) = t_hC;
for i=1:nreg
    theta(options.inF.indhD{i}) = t_hD{i};
end
theta(options.inF.indconst) = -.25;

%- observation
phi(options.inG.indr) = 0;

%--- Simulate time series of hidden states and observations
%disp('*** Simulation');
[y,x,x0,eta,e] = simulateNLSS(n_t,f_fname,g_fname,theta,phi,u,alpha,sigma,options);

% plot(x(options.inF.r,:))
% plot(y(2,:)'-e(2,:)')

% display time series of hidden states and observations
% displaySimulations(y,x,eta,e)
%  disp('--paused--')
%  pause

%-----------------------------------------------------------
%------------------- model inversion -----------------------
%--- Call inversion routine
disp('*** Hypothesis inversion');
options.graphic=0;
options.checkGrads = 0;
% theta

[posterior,out] = VBA_NLStateSpaceModel(y,u,f_fname,g_fname,dim,options);

% displayResults(posterior,out,y-e,x,x0,theta,phi,alpha,sigma)
% end