function [results] = demo_DCMdepth(n,nr)

if nargin==0
    n=2;
end

%-----------------------------------------------------------
%-------------- DCM model specification --------------------

% === Basic settings =======================================
f_fname = @f_DCMwHRFext ;  %
g_fname = @g_demo_extended;   % 

                              % 
TR = 2;                      % sampling period (in sec)
% n_t = round(15*60/TR);          % number of time samples
microDT = 2e-1;               % micro-time resolution (in sec)
homogeneous = 1;              % params of g(x) homogeneous accross regions
reduced_f = 1;                % fix some HRF params
lin = 1;                      % linearized variant of HRF Balloon model
stochastic = 0;               % flag for stochastic DCM inversion
alpha = Inf;                  % state noise precision
sigma = 1/5;                % measurement noise precision

% === Input ================================================

repets=30;
isi = zeros(1,5);
% factorial_u = [isi [1;0] isi [0;1] isi [1;1]];
factorial_u = [isi 1 0 0 0 0 isi .5 0 0 0 0];
factorial_s = [isi 0 0 1 1 0 isi 0 0 1 1 0];
u = [isi repmat(factorial_u,1,repets) isi isi];  
sampling = find([isi repmat(factorial_s,1,repets) isi isi]);  

nu = size(u,1);
n_t=size(u,2);


% === DCM structure ========================================
% invariant effective connectivity
A=zeros(n);
A(2:end,1:end-1)= eye(n-1);
nreg = size(A,1); 
% modulatory effects
B{1} = [] ;
% input-state coupling
C = zeros(n,1);
C(1,1)=1;
% gating (nonlinear) effects
for i=1:n
    D{i} = zeros(nreg);
end
% === Decoding scheme ========================================
hA = zeros(1,n);
hA(nr)=1;
nrep=size(hA,1);
hB = {zeros(nrep,nreg)}; 
hC = zeros(nrep,nu);
for i=1:n
    hD{i} = zeros(nrep,nreg);
end


% === Build options and dim structures =====================
options = prepare_fullDCM(A,B,C,D,TR,microDT,homogeneous,hA,hB,hC,hD);

options.priors = getPriors(nreg,n_t,options,reduced_f,stochastic);
options.microU = 0;
options.backwardLag = 8;
options.GnFigs = 0;
options.inF.linearized = lin;
options.DisplayWin=0;
options.verbose=0;
dim.n_theta = options.inF.ind5(end);
if options.extended
    dim.n_phi = options.inG.indr;
else
    dim.n_phi = options.inG.ind2(end);
end
dim.n = 5*nreg+nrep;
dim.n_t=n_t;


sources(1) = struct('out',1:n,'type',0); % BOLD signal (gaussian, dim=4)
sources(2) = struct('out',n+1,'type',1); % first binary response
options.sources=sources;


%% -----------------------------------------------------------
%----------- simulated times series specification ----------

%--- simulated evolution parameters: neuronal level
theta = zeros(dim.n_theta,1);
theta(options.inF.indA)=1;
theta(options.inF.indC)=1;
theta(options.inF.indhA)=1;
theta(options.inF.indself) = 0;
theta(options.inF.indconst) = -.1;
theta(options.inF.indhself) = log(1/options.inF.deltat) ;

phi = zeros(dim.n_phi,1);


[y,x,x0,eta,e] = simulateNLSS(n_t,f_fname,g_fname,theta,phi,u,alpha,sigma,options);
%    
% results.y=y;
% results.e=e;
% results.x=x;
% return


parfor i=1:20
    [y,x,x0,eta,e] = simulateNLSS(n_t,f_fname,g_fname,theta,phi,u,alpha,sigma,options);

    F_WO = zeros(1,n);
    F_W = zeros(1,n);
    
    options_n=options;
    options_n.isYout=zeros(size(y));
    options_n.isYout(end,:)=1;
    [posterior_WO,out_WO] = VBA_NLStateSpaceModel(y,u,f_fname,g_fname,dim,options_n);
    
    F_WO(1)=out_WO.F;
    priors_WO = options_n.priors;
    for nn=1:n-1
        subA = options_n.inF.indA(end-nn+1:end);
        priors_WO.SigmaTheta(subA,subA) = 0;
        F_WO(nn+1) = VB_SavageDickey(posterior_WO,options_n.priors,F_WO(1),out_WO.dim,priors_WO);
    end
%     priors_WO.SigmaTheta(options_n.inF.indA,options_n.inF.indA) = 0;
%     F_WO(end+1) = VB_SavageDickey(posterior_WO,options_n.priors,F_WO(1),out_WO.dim,priors_WO);

    options_n.isYout(end,:)=0;
    [posterior_W,out_W] = VBA_NLStateSpaceModel(y,u,f_fname,g_fname,dim,options_n);

    F_W(1)=out_W.F;
    priors_W = options_n.priors;
    for nn=1:n-1
        subA = options_n.inF.indA(end-nn+1:end);
        priors_W.SigmaTheta(subA,subA) = 0;
        F_W(nn+1) = VB_SavageDickey(posterior_W,options_n.priors,F_W(1),out_W.dim,priors_W);
    end
%     priors_W.SigmaTheta(options_n.inF.indA,options_n.inF.indA) = 0;
%     F_W(end+1) = VB_SavageDickey(posterior_W,options_n.priors,F_W(1),out_W.dim,priors_W);
%     
    results(i).n=n;
    results(i).nr=nr;
    results(i).theta=theta;
    results(i).muTheta_W=posterior_W.muTheta;
    results(i).SigmaTheta_W=posterior_W.SigmaTheta;
    results(i).muTheta_WO=posterior_WO.muTheta;
    results(i).SigmaTheta_WO=posterior_WO.SigmaTheta;
    results(i).noise=1/sigma;
    results(i).repet=repets;
    results(i).F_W=F_W;
    results(i).F_WO=F_WO;
end



end