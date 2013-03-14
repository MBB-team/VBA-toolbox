function [posterior,out,y] = demo_gonogo()
%%% ----------------------------------------------------------------------
%   Initialise the generative model, simulate data and inverse the model
%   
%   IN:     - b ; idyosincratic value of connectivity
%   OUT: classical structure from the DCM toolbox
%%% ----------------------------------------------------------------------

close all

%-----------------------------------------------------------
%-------------- DCM model specification --------------------

% === Basic settings =======================================
f_fname = @f_DCMwHRFext ;  %
g_fname = @g_HRF3ext;   % 

                              % 
TR = 2e0;                     % sampling period (in sec)
n_t = round(300/TR);          % number of time samples
microDT = 2e-1;               % micro-time resolution (in sec)
homogeneous = 1;              % params of g(x) homogeneous accross regions
reduced_f = 1;                % fix some HRF params
lin = 1;                      % linearized variant of HRF Balloon model
stochastic = 0;               % flag for stochastic DCM inversion
alpha = Inf;                  % state noise precision
sigma = 1/2;                  % measurement noise precision

% === Input ================================================
u = zeros(2,n_t);
% 
for(i=10:16:n_t)
    u(1,i:i+6) = rand(1)  ;
end
for(i=10:32:n_t)
    u(2,i+2:i+8) = 1 ;
end

sampler = zeros(1,n_t);
for(i=10:16:n_t)
    sampler(1,i+5) = 1  ;
end

nu = size(u,1);

% === DCM structure ========================================
% invariant effective connectivity
A = [0 0;
     0 0];
nreg = size(A,1);
% modulatory effects
B{1} = zeros(nreg,nreg);
B{2} = zeros(nreg,nreg);
% input-state coupling
C = zeros(nreg,2);
C(1:nreg,:) = eye(nreg,2);
% gating (nonlinear) effects
D{1} = zeros(nreg,nreg);
D{2} = zeros(nreg,nreg);
D{3} = zeros(nreg,nreg);
D{4} = zeros(nreg,nreg);
% === Decoding scheme ========================================
hA = zeros(1,nreg);
hA = [1 1];

nrep=size(hA,1);
hB = {zeros(nrep,nreg),zeros(nrep,nreg)}; 
hC = zeros(nrep,nu);
hD = {  zeros(nrep,nreg),...
        zeros(nrep,nreg),...
};

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

%% -----------------------------------------------------------
%----------- simulated times series specification ----------

%--- simulated evolution parameters: neuronal level
theta = zeros(dim.n_theta,1);
phi = zeros(dim.n_phi,1);

%- DCM
t_Aself = -.5;
t_A = [];
t_B{1} = [];
t_B{2} = [] ; 
t_C = [.8 .8];
t_D{1} = [];
t_D{2} = [];
t_D{3} = [];
t_D{4} = [];

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
t_hA = [.8 0];
t_hB{1} = [];
t_hB{2} = [];
t_hC = [];
t_hD{1} = [];
t_hD{2} = [];
t_hself = log(options.inF.deltat) ;

theta(options.inF.indhA) = t_hA;
theta(options.inF.indhself) = t_hself;
for i=1:nu
    theta(options.inF.indhB{i}) = t_hB{i};
end
theta(options.inF.indhC) = t_hC;
for i=1:nreg
    theta(options.inF.indhD{i}) = t_hD{i};
end

%- observation
phi(options.inG.indr) = 0;

%--- Simulate time series of hidden states and observations
%disp('*** Simulation');

[y,x,x0,eta,e] = simulateNLSS(n_t,f_fname,g_fname,theta,phi,u,alpha,sigma,options);



% display time series of hidden states and observations
% displaySimulations(y,x,eta,e)
%  disp('--paused--')
%  pause

%-----------------------------------------------------------
%------------------- model inversion -----------------------
%--- Call inversion routine
disp('*** Hypothesis inversion');
% A = [0 0 0 0;
%      0 0 0 0;
%      1 1 0 1;
%      1 1 1 0];
%  options = prepare_fullDCM(A,B,C,D,TR,microDT,homogeneous,hA,hB,hC,hD);

options.priors = getPriors(nreg,n_t,options,reduced_f,stochastic);
options.microU = 0;
options.backwardLag = 8;
options.GnFigs = 0;
options.inF.linearized = lin;
options.DisplayWin=1;
options.verbose=1;
dim.n_theta = options.inF.ind5(end);

options.isYout=zeros(size(y));
%options.isYout(3,:)=1-sampler;
%y(1:2,:) = randn(2,n_t);

options.graphic=1;
% options.TolFun=.00001;
options.checkGrads = 0;

[posterior,out] = VBA_NLStateSpaceModel(y,u,f_fname,g_fname,dim,options);

displayResults(posterior,out,y-e,x,x0,theta,phi,alpha,sigma)
end