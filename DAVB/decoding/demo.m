function [y,posterior,out] = demo()
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
n_t = round(550/TR);          % number of time samples
microDT = 2e-1;               % micro-time resolution (in sec)
homogeneous = 0;              % params of g(x) homogeneous accross regions
reduced_f = 0;                % fix some HRF params
lin = 1;                      % linearized variant of HRF Balloon model
stochastic = 0;               % flag for stochastic DCM inversion
alpha = Inf;                  % state noise precision
sigma = Inf;                  % measurement noise precision

% === Input ================================================
u = zeros(2,n_t);
% 
for(i=10:30:260)
    u(1,i:i+5) = 1  ;
    u(2,i:i+5) = .5 ;
    u(1,i+15:i+20) = .5  ;
    u(2,i+15:i+20) = 1 ;
end


nu = size(u,1);

% === DCM structure ========================================
% invariant effective connectivity
A = [0 0 0 0;
     0 0 0 0;
     1 0 0 0;
     0 1 0 0];
nreg = size(A,1);
% modulatory effects
B{1} = zeros(nreg,nreg);
B{2} = zeros(nreg,nreg);
% input-state coupling
C = [1 0; 
     0 1;
     0 0;
     0 0];
% gating (nonlinear) effects
D{1} = zeros(nreg,nreg);
D{2} = zeros(nreg,nreg);
D{3} = zeros(nreg,nreg);
D{4} = zeros(nreg,nreg);
% === Decoding scheme ========================================
hA = [0 0 1 1;
      0 0 1 1];
nrep=size(hA,1);
hB = {zeros(nrep,nreg),zeros(nrep,nreg)}; 
hC = zeros(nrep,nu);
hD = {  zeros(nrep,nreg),...
        zeros(nrep,nreg),...
        zeros(nrep,nreg),...
        zeros(nrep,nreg) ...
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
t_A = [-.3 -.3];
t_B{1} = [];
t_B{2} = [] ; 
t_C = [.4 .4];
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
t_hA = 2*[.8 -.8 -.8 .8];
t_hB{1} = [];
t_hB{2} = [];
t_hC = [];
t_hD{1} = [];
t_hD{2} = [];
t_hD{3} = [];
t_hD{4} = [];
t_hself = -.4 ;

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
phi(options.inG.indr) = log(40);

%--- Simulate time series of hidden states and observations
%disp('*** Simulation');

x = NaN;
while isweird(x)
    [y,x,x0,eta,e] = simulateNLSS(...
        n_t,f_fname,g_fname,theta,phi,u,alpha,sigma,options);
end
%[y,x,x0,eta,e] = simulateNLSS(n_t,f_fname,g_fname,theta,phi,u,alpha,sigma,options); 

% trick to get binomial response

y(5:6,:) = binomial_sample(y(5:6,:));
options.isYout = zeros(size(y));
options.isYout(5:6,:) = 1;
options.isYout(:,15:15:270) = 0;
%y(5:6,options.isYout==1)=NaN;
[sum(y(5,options.isYout(5,:)==0)) sum(y(6,options.isYout(6,:)==0))]

% display time series of hidden states and observations
displaySimulations(y,x,eta,e)
% disp('--paused--')
% pause

%y = repmat(mean(y,2),1,size(y,2)) +  (randn(size(y))'*chol(diag(var(y')')))' ;
%-----------------------------------------------------------
%------------------- model inversion -----------------------
%--- Call inversion routine
%disp('*** Hypothesis inversion');
options.graphic=1;
[posterior,out] = VBA_NLStateSpaceModel(y,u,f_fname,g_fname,dim,options);


end