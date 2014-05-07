function [posterior,out,nodes] = demo_susceptibility2()
%%% ----------------------------------------------------------------------
%   Initialise the generative model, simulate data and inverse the model
%   
%   IN:     feedback 
%%% ----------------------------------------------------------------------

% close all

%-----------------------------------------------------------
%-------------- DCM model specification --------------------

% === Basic settings =======================================
f_fname = @f_DCMwHRFext ;  %
g_fname = @g_demo_susceptibility;   % 

                              % 
TR = 1;                      % sampling period (in sec)
microDT = .1;               % micro-time resolution (in sec)
homogeneous = 1;              % params of g(x) homogeneous accross regions
reduced_f = 1;                % fix some HRF params
lin = 1;                      % linearized variant of HRF Balloon model
stochastic = 0;               % flag for stochastic DCM inversion
alpha = Inf;                  % state noise precision
sigma = [1/.05 1/.05];                % measurement noise precision

% === Input ================================================
nu = 3 ;
reps = 25;
u=zeros(nu,8*reps*(2^nu));

u_temp = [];
u_perm = full(spm_perm_mtx(nu))' ;
for i=1:reps
    u_temp = [u_temp u_perm(:,randperm(2^nu))];
end
u(:,1:8:end) = u_temp;
u(:,2:8:end) = u_temp;
u(:,3:8:end) = u_temp;

isYout=ones(5,size(u,2));
isYout(1:4,1:2:end) = 0;
% isYout(5,1:8:end) = 0;
% isYout(5,2:8:end) = 0;
isYout(5,3:8:end) = 0;

n_t=size(u,2);


% === DCM structure ========================================
% invariant effective connectivity

A = [0 0 0 0;
     1 0 0 0;
     0 0 0 0;
     0 0 0 0];

nreg = size(A,1);
% modulatory effects
B{1} = zeros(nreg,nreg);
B{2} = [0 0 0 0;
        0 0 0 0;
        0 1 0 0;
        0 0 0 0];
B{3} = [0 0 0 0;
        0 0 0 0;
        0 0 0 0;
        0 1 0 0];
% input-state coupling
C = [1 0 0;
     0 0 0;
     0 0 0;
     0 0 0];

% gating (nonlinear) effects
for i=1:nreg
D{i} = zeros(nreg,nreg);
end
% === Decoding scheme ========================================
hA = [0 0 1 1];
nrep=size(hA,1);
hB = {zeros(nrep,nreg),zeros(nrep,nreg)}; 
hC = zeros(nrep,nu);
hD = {  zeros(nrep,nreg),...
        zeros(nrep,nreg) ...
     };

% === Build options and dim structures =====================
options = prepare_fullDCM(A,B,C,D,TR,microDT,homogeneous,hA,hB,hC,hD);

sources(1) = struct('out',1:4,'type',0); % BOLD signal (gaussian, dim=4)
sources(2) = struct('out',5,'type',0); % first binary response
options.sources=sources;


options.priors = getPriors(nreg,n_t,options,reduced_f,stochastic);
options.microU = 0;
options.backwardLag = 8;
options.GnFigs = 0;
options.inF.linearized = lin;
options.DisplayWin=1;
options.verbose=0;


options.isYout = isYout ;
dim.n_theta = options.inF.ind5(end);
if options.extended
    dim.n_phi = options.inG.indr;
else
    dim.n_phi = options.inG.ind2(end);
end
dim.n = 5*nreg+nrep;
dim.n_t=n_t;




%% -----------------------------------------------------------
%----------- simulated times series specification ----------

%--- simulated evolution parameters: neuronal level
theta = zeros(dim.n_theta,1);
phi = zeros(dim.n_phi,1);

%- DCM
t_Aself = log(1);

t_A = [1];


    
t_B{1} = []; 
t_B{2} = 1; 
t_B{3} = 1; 

t_C = 1;
for i=1:nreg, t_D{i} = []; end


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
for i=1:nu; t_hB{i} = []; end
t_hC = [];
for i=1:nreg; t_hD{i} = []; end

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

theta(options.inF.indconst) = 0;


%- observation
% phi(options.inG.indr) = 100;

%--- Simulate time series of hidden states and observations
%disp('*** Simulation');
[y,x,x0,eta,e] = simulateNLSS(n_t,f_fname,g_fname,theta,phi,u,alpha,sigma,options);
% figure ;
% plot(y(1:4,:)'-e(1:4,:)')
% posterior=[];
% out=[];
% return
% plot(y');
%-----------------------------------------------------------
%------------------- model inversion -----------------------
%--- Call inversion routine
% disp('*** Hypothesis inversion');


options = prepare_fullDCM(A,B,C,D,TR,microDT,homogeneous,hA,hB,hC,hD);
options.sources=sources;

options.priors = getPriors(nreg,n_t,options,reduced_f,stochastic);
options.microU = 0;
options.backwardLag = 8;
options.GnFigs = 0;
options.inF.linearized = lin;
options.DisplayWin=1;
options.verbose=0;
options.graphic=0;
options.checkGrads = 0;
dim.n_theta = options.inF.ind5(end);
if options.extended
    dim.n_phi = options.inG.indr;
else
    dim.n_phi = options.inG.ind2(end);
end
dim.n = 5*nreg+nrep;
dim.n_t=n_t;

options.isYout = isYout ;


[posterior,out] = VBA_NLStateSpaceModel(y,u,f_fname,g_fname,dim,options);

%%
nodes= grapher_getDcmTemplate([0 150; 0 -20; -150 -150; 150 -150],{'1','2','3','4'});


end