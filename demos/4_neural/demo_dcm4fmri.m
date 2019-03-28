% Demo for (stochastic) DCM for fMRI.
% This demo inverts the DCM for fMRI model, which contains the ballon model
% as a generalized observation function (not affected by stochastic
% innovations).
% NB: currently, there is a ratio 1e2 between the expected standard
% deviations of neural versus hemodynamic (states) noise. This means the
% emphasis is on recovering neural noise, which partly drives the neural
% dynamics and is passed through the HRF dynamical model.

close all
clear variables

%-----------------------------------------------------------
%-------------- DCM model specification --------------------

%--- Basic settings
f_fname = @f_DCMwHRF;
g_fname = @g_HRF3;
TR = 2e0;                     % sampling period (in sec)
n_t = round(1e2/TR);          % number of time samples
dtU = round(2/TR)+1;          % input-on time interval
t0U = round(10/TR)+1;                  
microDT = 2e-1;               % micro-time resolution (in sec)
homogeneous = 0;              % params of g(x) homogeneous accross regions
reduced_f = 0;                % fix some HRF params
lin = 1;                      % linearized variant of HRF Balloon model
stochastic = 0;               % flag for stochastic DCM inversion
alpha = Inf;%1e2/TR;               % state noise precision
sigma = 1e0;                  % measurement noise precision
nconfounds = 0;               % # of basis sets in the confounds matrix

%--- Input
u       = zeros(2,n_t);
u(1,t0U:t0U+dtU) = 1;
u(1,4*t0U:4*t0U+dtU) = 1;
u(2,4*t0U:4*t0U+5*dtU) = 1;
nu = size(u,1);

%--- DCM structure
% invariant effective connectivity
A = [0 1 1
     1 0 1
     0 1 0];
nreg = size(A,1);
% modulatory effects
B{1} = zeros(nreg,nreg);
B{2} = [0 0 0
        1 0 0
        0 0 0];
% input-state coupling
C = [1 0
     0 0
     0 0];
% gating (nonlinear) effects
D{1} = [0 0 0
        0 0 0
        0 0 0];
D{2} = zeros(nreg,nreg);
D{3} = zeros(nreg,nreg);

%--- Build options and dim structures
options = prepare_fullDCM(A,B,C,D,TR,microDT,homogeneous);
options.priors = getPriors(nreg,n_t,options,reduced_f,stochastic);
options.microU = 0;
options.backwardLag = ceil(16/TR);  % 16 secs effective backward lag
options.inF.linearized = lin;
dim.n_theta = options.inF.ind5(end);
if options.inG.homogeneous
    dim.n_phi = 2;
else
    dim.n_phi = 2*nreg;
end  
dim.n = 5*nreg;


%-----------------------------------------------------------
%----------- simulated times series specification ----------

%--- simulated evolution parameters: neuronal level
t_A = exp([ -0.5
            -0.5
            -0.5
            -2.5
            -1.5 ]);
t_Aself = -0;
t_B{1} = [];
t_B{2} = exp([ -0.5 ]);
t_C = exp([ +0.1 ]);
t_D{1} = 0;%-exp([ -1 ]);
t_D{2} = [];
t_D{3} = [];
theta = zeros(dim.n_theta,1);
phi = zeros(dim.n_phi,1);
theta(options.inF.indA) = t_A;
for i=1:nu
    theta(options.inF.indB{i}) = t_B{i};
end
theta(options.inF.indC) = t_C;
for i=1:nreg
    theta(options.inF.indD{i}) = t_D{i};
end


%--- Simulate time series of hidden states and observations
x = NaN;
while VBA_isWeird (x)
    [y,x,x0,eta,e] = VBA_simulate (n_t,f_fname,g_fname,theta,phi,u,alpha,sigma,options);
end

% add in confounds
if nconfounds > 0
    [X0] = get_U_basis(n_t*TR,TR,nconfounds,'Fourier')';
    P0 = randn(nconfounds*nreg,1);
    phi = [phi;P0];
    P0 = reshape(P0,nconfounds,nreg);
    y = y + [X0*P0]';
%     figure,plot([X0*P0])
    [u,options,dim] = addConfounds2dcm(X0,u,options,dim);
end



% display time series of hidden states and observations
displaySimulations(y,x,eta,e);
% disp('--paused--')
% pause

% options.isYout = zeros(size(y));
% for t=1:n_t
%     options.priors.iQy{t} = eye(3);
%     if t > 20 && t < 50
% %        options.priors.iQy{t} = 0.*options.priors.iQy{t};
%         options.isYout(:,t) = 1;
%     end
% end


%-----------------------------------------------------------
%------------------- model inversion -----------------------
%--- Call inversion routine
% options.checkGrads = 1;
% options.gradF = 1;
% options.MaxIter = 2;
[posterior,out] = VBA_NLStateSpaceModel(y,u,f_fname,g_fname,dim,options);

%--- Display results
displayResults(posterior,out,y-e,x,x0,theta,phi,alpha,sigma);

%--- Make predictions
try
    options = out.options;
    [xs,ys,xhat,vx,yhat,vy] = VBA_comparePredictions(n_t,theta,phi,u,alpha,sigma,options,posterior,dim);
catch
    disp('------!!Unable to form predictions!!------')
end


