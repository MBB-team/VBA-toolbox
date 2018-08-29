% Demo for DAVB Henon state-space model inversion...

clear variables
close all

% Choose basic settings for simulations
n_t = 2e2;
f_fname = @f_Henon;
g_fname = @g_sigmoid;
u       = [];
alpha   = 1e8;
sigma   = 1e4;
theta   = [1.4;0.3];
phi     = [];
inF.deltat = 1;
inG.scale = 50;
inG.slope = 0.1;

% Build options structure for temporal integration of SDE
options.inF     = inF;
options.inG     = struct;
% options.u0      = 0*ones(2,1);  % initial condition: input value


% Build priors for model inversion

priors.muX0 = 0*ones(2,1);
priors.SigmaX0 = 1e0*eye(2);

priors.muTheta = 0.*ones(2,1);
priors.SigmaTheta = 1e0*eye(2);

% priors.muPhi = zeros(2,1);
% priors.SigmaPhi = 1e4*speye(2);

priors.a_alpha = 1;
priors.b_alpha = 1;

priors.a_sigma = 1e0;
priors.b_sigma = 1e0;


options.priors      = priors;
options.backwardLag = 4;
dim.n_theta         = 2;
dim.n_phi           = 0;
dim.n               = 2;



% Build time series of hidden states and observations
ok = false;
while ~ ok
    try
        [y,x,x0,eta,e] = VBA_simulate (n_t,f_fname,g_fname,theta,phi,u,alpha,sigma,options);
        ok = true;
    catch
        fprintf('Simulation produced weird values, let''s try again.\n');
    end
end


% display time series of hidden states and observations
displaySimulations(y,x,eta,e);
% disp('--paused--')
% pause

% options.checkGrads = 1;

% Call inversion routine
[posterior,out] = VBA_NLStateSpaceModel(y,u,f_fname,g_fname,dim,options);


%------------ Display results ------------------%
displayResults(posterior,out,y,x,x0,theta,phi,alpha,sigma);

% Make predictions
try
    options = out.options;
    [xs,ys,xhat,vx,yhat,vy] = VBA_comparePredictions(...
        n_t,theta,phi,u,alpha,sigma,options,posterior,dim);
catch
    disp('------!!Unable to form predictions!!------')
end



