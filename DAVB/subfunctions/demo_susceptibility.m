function [posterior,out,t_A] = demo_susceptibility(nu,reps,recurs)
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

u=zeros(nu,4*reps*(2^nu));

u_temp = [];
u_perm = full(spm_perm_mtx(nu))' ;
for i=1:reps
    u_temp = [u_temp u_perm(:,randperm(2^nu))];
end
u(:,1:4:end) = u_temp;

n_t=size(u,2);


% === DCM structure ========================================
% invariant effective connectivity


if recurs
    A_temp = [ zeros(1,nu) ; eye(nu-1,nu) ];
    A = [zeros(nu)  A_temp; eye(nu) zeros(nu)];
else
    A = [zeros(nu)  zeros(nu) ; ones(nu) zeros(nu)];
end

nreg = size(A,1);
% modulatory effects
for i=1:nu
    B{i} = zeros(nreg,nreg);
end
% input-state coupling
C = [eye(nu) ; zeros(nu)];

% gating (nonlinear) effects
for i=1:nreg
D{i} = zeros(nreg,nreg);
end
% === Decoding scheme ========================================
hA = [zeros(1,nu) ones(1,nu)];
nrep=size(hA,1);
hB = {zeros(nrep,nreg),zeros(nrep,nreg)}; 
hC = zeros(nrep,nu);
hD = {  zeros(nrep,nreg),...
        zeros(nrep,nreg) ...
     };

% === Build options and dim structures =====================
options = prepare_fullDCM(A,B,C,D,TR,microDT,homogeneous,hA,hB,hC,hD);

sources(1) = struct('out',1:(2*nu),'type',0); % BOLD signal (gaussian, dim=4)
sources(2) = struct('out',(2*nu)+1,'type',0); % first binary response
options.sources=sources;


options.priors = getPriors(nreg,n_t,options,reduced_f,stochastic);
options.microU = 0;
options.backwardLag = 8;
options.GnFigs = 0;
options.inF.linearized = lin;
options.DisplayWin=1;
options.verbose=1;
options.isYout=ones(nu^2+1,size(u,2));
options.isYout(1:end-1,1:2:end) = 0;
options.isYout(end,1:4:end) = 0;

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

if recurs
    t_A = randi(3,1,2*nu-1);
else
    t_A = randi(3,1,nu^2 );
end

    
for i=1:nu, t_B{i} = []; end
t_C = ones(1,nu);
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

t_hA = (4/nu)*ones(1,nu); 
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
% plot(y(end,:)-e(end,:))

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
options.verbose=1;
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


options.isYout=ones(nu^2+1,size(u,2));
options.isYout(1:end-1,1:2:end) = 0;
options.isYout(end,1:4:end) = 0;

[posterior,out] = VBA_NLStateSpaceModel(y,u,f_fname,g_fname,dim,options);



end