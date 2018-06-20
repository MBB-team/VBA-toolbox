% This demonstrates the smoothing properties of the Kalman lagged filter.
% Hidden states follow a triangular wave, whose observation is perturbed
% with white noise. Critically, we render the inversion scheme blind during
% half a period of the oscillation. We then invert the model (under AR
% priors on hidden states), with and without large backward lag.

clear variables
close all

nt = 10;
eta = 0.*randn(1,3*nt);
x = 1:nt;
x = [x,fliplr(x),x] + eta;

e = randn(1,3*nt);
y = x+e;

theta = [];
phi = [];
alpha = 1/var(eta);
sigma = 1/var(e);

displaySimulations(y,x,eta,e);


f_fname = @f_AR;
g_fname = @g_Id;
dim.n_theta = size(theta,1);
dim.n_phi = size(phi,1);
dim.n = size(x,1);

options.isYout = zeros(1,3*nt);
options.isYout(nt:2*nt+1) = 1;
options.MaxIterInit = 0;
options.priors.a_alpha = 1;
options.priors.b_alpha = 1;

% VB-Kalman-filter
[p1,o1] = VBA_NLStateSpaceModel(y,[],f_fname,g_fname,dim,options);

% VB-Kalman-smoother (lag = size of blind window)
options.backwardLag = nt+1;
[p2,o2] = VBA_NLStateSpaceModel(y,[],f_fname,g_fname,dim,options);
set(gcf,'name',['Kalman lag = ',num2str(o2.options.backwardLag)])

VBA_ReDisplay(p1,o1,1);
set(gcf,'name',['Kalman lag = ',num2str(o1.options.backwardLag)])