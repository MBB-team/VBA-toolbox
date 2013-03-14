function [posterior,out] = demo_extended(model)
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
g_fname = @g_demo_extended;   % 

                              % 
TR = 2e0;                     % sampling period (in sec)
n_t = round(175/TR);          % number of time samples
microDT = 1e-1;               % micro-time resolution (in sec)
homogeneous = 1;              % params of g(x) homogeneous accross regions
reduced_f = 0;                % fix some HRF params
lin = 1;                      % linearized variant of HRF Balloon model
stochastic = 0;               % flag for stochastic DCM inversion
alpha = Inf;                  % state noise precision
sigma = Inf;                  % measurement noise precision

% === Input ================================================
% u = zeros(2,n_t);
% % 
% for(i=10:32:n_t)
%     u(1,i:i+8) = 1  ;
%     u(2,i:i+8) = .5 ;
%     u(1,i+16:i+24) = .5  ;
%     u(2,i+16:i+24) = 1 ;
% end
% 
% 
% nu = size(u,1);

 u = zeros(2,n_t);
% orthogonal modulation
cpt_0=0;
cpt_1=0;
for(i=10:20:n_t)
    if mod(cpt_0,2) <1
        u(1,i:i+8) = 1  ;
    else
        u(1,i:i+8) = .8  ;
    end
    if mod(cpt_1,4) < 2
        u(2,i:i+8) = 1  ;
    else
        u(2,i:i+8) = .8  ;
    end
    cpt_0=cpt_0+1;
    cpt_1=cpt_1+1;
    
end


nu = size(u,1);


% === DCM structure ========================================
% invariant effective connectivity
switch model
    case 'race'
A = [0 0 0 0;
     0 0 0 0;
     1 0 0 0;
     0 1 0 0];
    case 'forward' 
A = [0 0 0 0;
     0 0 0 0;
     1 1 0 0;
     1 1 0 0];
    case 'lateral' 
A = [0 0 0 0;
     0 0 0 0;
     1 0 0 1;
     0 1 1 0];
    otherwise
        error('*** model should be: race, forward or lateral\n') ;
end
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
hA = zeros(2,nreg);
hA(:,3:4) = eye(2,2);


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

sources(1) = struct('out',1:4,'type',0); % BOLD signal (gaussian, dim=4)
sources(2) = struct('out',5:6,'type',1); % first binary response
options.sources=sources;


%% -----------------------------------------------------------
%----------- simulated times series specification ----------

%--- simulated evolution parameters: neuronal level
theta = zeros(dim.n_theta,1);
phi = zeros(dim.n_phi,1);

%- DCM
t_Aself = .3;
switch model
    case 'race'
    t_A = [.6 .6];
    case 'forward' 
	t_A = [.5 -.5 -.5 .5];
    case 'lateral'
    t_A = [.7 .7 -1 -1];
end
t_B{1} = [];
t_B{2} = [] ; 
t_C = [1 1];
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

t_hA = [1 1];

t_hB{1} = [];
t_hB{2} = [];
t_hC = [];
t_hD{1} = [];
t_hD{2} = [];
t_hD{3} = [];
t_hD{4} = [];
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
A = [0 0 0 0;
     0 0 0 0;
     1 1 0 1;
     1 1 1 0];
options = prepare_fullDCM(A,B,C,D,TR,microDT,homogeneous,hA,hB,hC,hD);
options.priors = getPriors(nreg,n_t,options,reduced_f,stochastic);
options.microU = 0;
options.backwardLag = 8;
options.GnFigs = 0;
options.inF.linearized = lin;
options.DisplayWin=0;
options.verbose=1;
dim.n_theta = options.inF.ind5(end);
options.sources=sources;
% options.isYout=zeros(size(y));
%options.isYout(1,:)=1;
%y(1:2,:) = randn(2,n_t);

options.graphic=1;
options.checkGrads = 0;

[posterior,out] = VBA_NLStateSpaceModel(y,u,f_fname,g_fname,dim,options);

% displayResults(posterior,out,y-e,x,x0,theta,phi,alpha,sigma)
end