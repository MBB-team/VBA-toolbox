function [posterior,out,theta] = demo_negfeedback(connectivity_model, encoding_region, noise,nRepetition)
%%% ----------------------------------------------------------------------
%   Initialise the generative model, simulate data and inverse the model
%   
%   IN:     feedback 
%%% ----------------------------------------------------------------------

if nargin==0
    connectivity_model = 'feedback_modulation' ;
    encoding_region    = 'first' ;
    noise = 0;
    nRepetition = 10;
end

%% #########################################################
%     DCM model specification 
%  #########################################################

switch connectivity_model
    case 'forward_modulation'
       As = [ 0 -1 ];
       Bs = [1.5 0 ] ;
    case 'feedback_modulation'
       As = [1   0 ];
       Bs = [0 -1.5] ;
    case 'no_modulation'
       As = [1  -1 ];
       Bs = [0   0 ] ;
    otherwise 
       error('connectivity_model: ''forward_modulation'', ''feedback_modulation'', ''no_modulation''');
end

switch encoding_region
    case 'first'
       hAs = [2 0];  
    case 'second'
       hAs = [0 2];    
    otherwise 
       error('encoding_region: ''first'', ''second''');
end
% __________________________________________________________
% === Input ================================================

u_condOn  = [1 0 0 0 0 0 0 0;
             1 1 1 1 1 1 1 1];  
         
u_condOff = [1 0 0 0 0 0 0 0;
             0 0 0 0 0 0 0 0];
         
isi = zeros(2,8);

u = [isi];
for iRepet = 1:nRepetition
    jitter = zeros(2,randi(3)-1);
    u = [u u_condOn jitter isi] ;
    jitter = zeros(2,randi(3)-1);
    u = [u u_condOff jitter isi] ;
end

u = [u isi];
[nu, n_t] = size(u);

% __________________________________________________________
% === DCM structure ========================================

% invariant effective connectivity
A = [0 1;
     1 1];
nreg = size(A,1);
% modulatory effects
B{2} = [0 1;
        1 0];
% input-state coupling
C = [1 0;
     0 0];

% __________________________________________________________
% === Decoding scheme ======================================

hA = [1 1];

% __________________________________________________________
% === Basic settings =======================================

f_fname = @f_DCMwHRFext ; 
g_fname = @g_DCMwHRFext ; 
                              
TR = 1;                       % sampling period (in sec)
microDT = .200;                % micro-time resolution (in sec)
homogeneous = 1;              % params of g(x) homogeneous accross regions
reduced_f = 1;                % fix some HRF params
lin = 1;                      % linearized variant of HRF Balloon model
stochastic = 0;               % flag for stochastic DCM inversion
alpha = Inf;                  % state noise precision
sigma = [1/noise];              % measurement noise precision

% __________________________________________________________ 
% === Build options and dim structures =====================

%- specify distribution of observations
sources(1) = struct('out',1:2,'type',0);  % two BOLD timeseries (gaussian) 
sources(2) = struct('out',3,  'type',1);  % and a motor response(gaussian)

%- prepare the DCM structure
options = prepare_fullDCM(A,B,C,{},TR,microDT,homogeneous,hA,{[],[]},{},{},sources);

options.priors = getPriors(nreg,n_t,options,reduced_f,stochastic);
options.backwardLag = 8;
options.inF.linearized = lin;

options.dim.n_t = n_t;
dim=options.dim;

%- indicate the datapoints of interest
TRidx   = find(mod(1:n_t,2/TR)) ; % keep one BOLD image per TR
respIdx = find(u(1,:))          ; % keep 4 responses per trial
respIdx = sort([respIdx, respIdx+1, respIdx+2, respIdx+3]) ;   

options.isYout = ones(3,n_t);
options.isYout(1:2,TRidx)   = 0;
options.isYout(  3,respIdx) = 0 ;


%% #########################################################
%     Simulation 
%  #########################################################

% __________________________________________________________
% === simulated times series specification =================

theta = zeros(dim.n_theta,1);
phi = zeros(dim.n_phi,1);

%- DCM
theta(options.inF.indA)     = [As .4]; % intrinsic connectivity
theta(options.inF.indself)  = log(1);  
theta(options.inF.indB{2})  = Bs;      % modulatory influence (PPI)
theta(options.inF.indC)     = 1;       % direct inputs

%- decoding 
theta(options.inF.indhA)    = hAs;
theta(options.inF.indhself) = log(1/options.inF.deltat); %/options.inF.deltat

% __________________________________________________________
% === Simulate time series of hidden states and observations 

disp('*** Simulation');
[y,x,x0,eta,e,u] = VBA_simulate (n_t,f_fname,g_fname,theta,phi,u,alpha,sigma,options);

[hf] = displaySimulations(y,x,eta,e)

%% #########################################################
%     Inversion 
%  #########################################################

%-----------------------------------------------------------
%------------------- model inversion -----------------------
disp('*** model inversion');

% specify priors
options.priors = getPriors(nreg,n_t,options,reduced_f,stochastic);

% fix unused paramters
idd = [options.inF.indA, options.inF.indB{2}, options.inF.indhA];
options.priors.SigmaTheta(idd,idd) = options.priors.SigmaTheta(idd,idd) .* diag(theta(idd)~=0);

[posterior,out] = VBA_NLStateSpaceModel(y,u,f_fname,g_fname,dim,options);


end