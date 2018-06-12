% Demo for dummy HRF model inversion.
% This demo simulates an HRF using a Fourier basis set and parameters that
% have been estimated on real data.
% The demo then inverts a model based on a mixture of two Gamma functions,
% estimating parameters and deriving the model evidence.


close all
clear variables


%---- HRF estimation  ----%

% Choose basic settings for simulations
n_t = 1;
p = 60;     % the fMRI time series length (in scans)
f_fname = [];
g_fname = @g_Fourier;
theta = [];
phi = [ 0.0126
       -0.0157
        0.0012
       -0.0020
        0.0080
        0.0113
        0.0177
        0.0137
        0.0130
        0.0088
        0.0087
        0.0060
        0.0064
        0.0044
        0.0050
        0.0035 ];

alpha = [];
sigma = 1e5;

% Build options structure for temporal integration of SDE
options.delta_t = 0;         % integration time step (Euler method)
options.t0 = 0;                 % initial condition: time value
options.x0 = 0;        % initial condition: states value
options.inG.grid = [0.01:p-1+0.01]';
inG.p = p;

% Build time series of fMRI observations 
[gx,dG_dX,dG_dPhi] = feval(g_fname,[],phi,[],inG);
y = gx + sqrt(sigma.^-1)*randn(size(gx));

    
% display time series of hidden states and observations
figure,
plot(y','ro'),title('y : as build using Fourier basis set')
hold on;
plot(gx')
% disp('--paused--')
% pause

% Build priors for model inversion

% priors.muPhi =0*ones(size(phi)-[5,0]);
priors.muPhi = [log(5);log(1);log(7);log(1);1;-0.2;0];
priors.SigmaPhi = 1e-1*eye(size(priors.muPhi,1));

y0 = (1./sqrt(sigma))*randn(size(y));
priors.a_sigma = numel(y0);
priors.b_sigma = trace(y0'*y0);


options.priors = priors;
g_fname = @g_DoubleGamma;

dim.n_theta = 0;
dim.n_phi = 7;
dim.n=0;
u = 0;

% options.checkGrads = 1;

% Call inversion routine
[posterior,out] = VBA_NLStateSpaceModel(y,u,[],g_fname,dim,options);


