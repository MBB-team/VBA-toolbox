function [results] = demo_sensitivity()

%-----------------------------------------------------------
%-------------- DCM model specification --------------------

% === Basic settings =======================================
f_fname = @f_DCMwHRFext ;  %
g_fname = @g_demo_extended;   % 

                              % 
TR = 2;                      % sampling period (in sec)
% n_t = round(15*60/TR);          % number of time samples
microDT = 2e-1;               % micro-time resolution (in sec)
homogeneous = 0;              % params of g(x) homogeneous accross regions
reduced_f = 0;                % fix some HRF params
lin = 1;                      % linearized variant of HRF Balloon model
stochastic = 0;               % flag for stochastic DCM inversion
alpha = Inf;                  % state noise precision
sigma = 1/5;                % measurement noise precision

% === Input ================================================

repets=10;
isi = zeros(2,9);
factorial_u = [isi [1;0] isi [0;1] isi [1;1]];
u = [isi repmat(factorial_u,1,repets) isi isi];  
sampling = find(sum(u)>0) ;

nu = size(u,1);
n_t=size(u,2);


% === DCM structure ========================================
% invariant effective connectivity
A = [0 1;
     1 0];

nreg = size(A,1); 
% modulatory effects
B{1} = [0 1; 
        1 0]; 
B{2} = zeros(nreg);
% input-state coupling
C = [1 0;
     0 1];
% gating (nonlinear) effects
D{1} = zeros(nreg);
D{2} = zeros(nreg);
% === Decoding scheme ========================================
hA = [1 1];
nrep=size(hA,1);
hB = {zeros(nrep,nreg),zeros(nrep,nreg)}; 
hC = zeros(nrep,nu);
hD = {  zeros(nrep,nreg),...
        zeros(nrep,nreg) ...
     };

% === Build options and dim structures =====================
options = prepare_fullDCM(A,B,C,D,TR,microDT,homogeneous,hA,hB,hC,hD);

options.priors = getPriors(nreg,n_t,options,reduced_f,stochastic);
options.microU = 0;
options.backwardLag = 8;
options.GnFigs = 0;
options.inF.linearized = lin;
options.DisplayWin=0;
options.verbose=1;
dim.n_theta = options.inF.ind5(end);
if options.extended
    dim.n_phi = options.inG.indr;
else
    dim.n_phi = options.inG.ind2(end);
end
dim.n = 5*nreg+nrep;
dim.n_t=n_t;


sources(1) = struct('out',1:2,'type',0); % BOLD signal (gaussian, dim=4)
sources(2) = struct('out',3,'type',1); % first binary response
options.sources=sources;


%% -----------------------------------------------------------
%----------- simulated times series specification ----------

%--- simulated evolution parameters: neuronal level
theta = zeros(dim.n_theta,1);
theta_DCM = [options.inF.indA options.inF.indB{1} options.inF.indC];
theta_ext = [options.inF.indhA  options.inF.indhB{2}];


phi = zeros(dim.n_phi,1);
t_Aself = 0+(rand(1)-.5);
t_hself = log(1/options.inF.deltat) ;
theta(options.inF.indself) = t_Aself;
theta(options.inF.indhself) = t_hself;

effect = [0 0 0];
options.verbose=0;
while ~all(effect >= [0 0 .6] & effect <= [5 5 1])
    theta(theta_DCM) = .4*ones(1,length(theta_DCM),1).*(randi(3,1,length(theta_DCM))-2);
    theta(theta_ext) = 1*ones(1,length(theta_ext),1).*(randi(3,1,length(theta_ext))-2);
    [y,x,x0,eta,e] = simulateNLSS(50,f_fname,g_fname,theta,phi,u,alpha,sigma,options);
    effect = (max(y'-e')-min(y'-e'));
end

parfor n=1:50
    [y,x,x0,eta,e] = simulateNLSS(n_t,f_fname,g_fname,theta,phi,u,alpha,sigma,options);


    options_n=options;
    options_n.isYout=zeros(size(y));
    options_n.isYout(3,:)=1;
    [posterior_WO,out_WO] = VBA_NLStateSpaceModel(y,u,f_fname,g_fname,dim,options_n);

    options_n.isYout(3,sampling)=0;
    [posterior_W,out_W] = VBA_NLStateSpaceModel(y,u,f_fname,g_fname,dim,options_n);

    results(n).theta=theta;
    results(n).effect=effect;
    results(n).muTheta_W=posterior_W.muTheta;
    results(n).SigmaTheta_W=posterior_W.SigmaTheta;
    results(n).muTheta_WO=posterior_WO.muTheta;
    results(n).SigmaTheta_WO=posterior_WO.SigmaTheta;
    results(n).noise=1/sigma;
    results(n).repet=repets;
end



end