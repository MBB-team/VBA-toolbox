% Demo for 2D neural field model.
% This demo inverts a linear 2D neural field model.
% NB: in this instance, the neural field is deterministic.

clear variables
close all


%---- Choose basic settings for simulations ----%

% dimensions
n_t             = 2e2;                      % number of time points
ns              = 5;                        % size of the square grid
deltat          = 5e-2;                     % time discretization
f_fname         = @f_2DneuralField;         % Evolution function
g_fname         = @g_Id;                    % Observation function

% Sinusoidal input to the approximate center of the neural field
u               = zeros(ns^2,n_t);
ohmega          = 0.91;
centre          = 0.5*ns*(ns-1);
u(centre-1:centre+1,2)     = 10;
u(centre-1:centre+1,:)     = repmat(1*sin(ohmega*deltat*[1:n_t]),3,1);

% parameters of the simulation
alpha           = Inf;%1e6;
sigma           = 1e3;%Inf;
theta           = [0.15;1;1];
phi             = [1;1;0.1];


% Build options structure for temporal integration of SDE
I = speye(ns,ns);
E = sparse(2:ns,1:ns-1,1,ns,ns);
D = E+E'-2*I;
inF.L           = kron(D,I)+kron(I,D);
inF.deltat      = deltat;
inG.ind         = 1:ns^2*2; % this makes states observable (not temporal derivatives)
options.inF     = inF;
options.inG     = inG;
% options.u0      = 0*ones(ns^2,1);



% Build priors for model inversion
priors.muX0 = 0*ones(ns^2*2,1);
priors.SigmaX0 = 0*speye(ns^2*2);%1e-1*speye(ns^2*2);
priors.muTheta = [.5;.5;.5];
priors.SigmaTheta = 1e0*speye(size(theta,1));
priors.muPhi = [0.1;0.9;0];
priors.SigmaPhi = 1e0*speye(size(phi,1));
priors.SigmaPhi(2,2) = 1e-4;
priors.a_alpha = Inf;%1e3;
priors.b_alpha = 0;%1e1;
priors.a_sigma = 1e3;
priors.b_sigma = 1e1;

% Build options and dim structures for model inversion
options.priors      = priors;
dim.n_theta         = size(theta,1);
dim.n_phi           = size(phi,1);
dim.n               = (ns^2)*2;



% Build time series of hidden states and observations
[y,x,x0,eta,e] = VBA_simulate (n_t,f_fname,g_fname,theta,phi,u,alpha,sigma,options,priors.muX0);

% display time series of hidden states and observations
displaySimulations(y,x,eta,e);
% disp('--paused--')
% pause


% Call inversion routine
% [posterior,out] = VBA_onlineWrapper(y,u,f_fname,g_fname,dim,options);
[posterior,out] = VBA_NLStateSpaceModel(y,u,f_fname,g_fname,dim,options);


% Display results
displayResults(posterior,out,y,x,x0,theta,phi,alpha,sigma);

% Make predictions
try
    options = out.options;
    [xs,ys,xhat,vx,yhat,vy] = VBA_comparePredictions(...
        n_t,theta,phi,zeros(size(u)),alpha,sigma,options,posterior,dim);
catch
    disp('------!!Unable to form predictions!!------')
end

