% this demo inverts a dynamical system that is sampled on an irregular grid

clear all
close all

nt = 2e2;
f_fname = @f_2dwu;
g_fname = @g_Id;

in.A = [-4 -16
        4   -4];
in.dt = 1e-1;
x0 = [1;1];
theta = [1;2];
phi = [];
alpha = Inf;
sigma = 1;
u =randn(2,nt);

% Build options and dim structures for model inversion
options.inF = in;
options.inG = in;
dim.n_theta = 2;
dim.n_phi = 0;
dim.n = 2;
dim.p = 2;

% Build time series of hidden states and observations
[y,x,x0,eta,e] = simulateNLSS(nt,f_fname,g_fname,theta,phi,u,alpha,sigma,options,x0);
displaySimulations(y,x,eta,e)

% Invert deterministic model
options.priors.SigmaTheta = 1e2*eye(2);
[posterior,out] = VBA_NLStateSpaceModel(y,u,f_fname,g_fname,dim,options);
displayResults(posterior,out,y-e,x,x0,theta,phi,alpha,sigma);


% Now resample data on an irregular grid by 'removing' data
in = zeros(1,nt);
in(2.^[1:floor(log(nt)./log(2))]) = 1; % exponential sampling grid
isout = repmat(1-in,size(y,1),1); % point sto remove
y(isout==1) = 0;
options.isYout = isout;
[posterior,out] = VBA_NLStateSpaceModel(y,u,f_fname,g_fname,dim,options);
displayResults(posterior,out,y-e,x,x0,theta,phi,alpha,sigma);