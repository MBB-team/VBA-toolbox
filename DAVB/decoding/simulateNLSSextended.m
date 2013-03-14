function [y,x,r,x0,r0,eta,e] = simulateNLSSextended(n_t,f_fname,g_fname,h_fname,theta,phi,psi,u,alpha,sigma,options,x0,r0)
% samples times series from sDCM generative model
% [y,x,r,x0,dTime,eta,eta0] =
% simulateNLSS(n_t,f_fname,g_fname,h_fname,theta,phi,psi,u,alpha,sigma,options,x0)
%
% This function creates the time series of hidden-states and measurements
% under the following nonlinear state-space model:
%   x_t = f(x_t-1,Theta,u_t) + f_t
%   y_t = g(x_t,Phi,u_t) + e_t
%   r_t = h(x_t,Psi,u_t)

% where f, g, and h are the evolution, observation, and decoding functions respectively.
% IN:
%   - n_t: the number of time bins for the time series of hidden-states and
%   observations, i.e. the time indices satisfy: 1<= t < n_t
%   - f_fname/g_fname/hname: evolution/observation/decoding function names.
%   - theta/phi/psi: evolution/observation/decoding parameters values.
%   - u: the mxt input matrix
%   - alpha: precision of the stochastic innovations
%   - sigma: precision of the measurement error
%   - options: structure variable containing the following fields:
%       .inF
%       .inG
%       .inH
%   - x0: the initial conditions
% OUT:
%   - y: the pxt (noisy) measurement time series
%   - x: the nxt (noisy) hidden-states time series
%   - r: the nxt predictions of behavioral responses (Poisson parameter
%        time series)
%   - x0: the nx1 initial conditions
%   - eta: the nxt stochastic innovations time series
%   - e: the pxt measurement errors (e:=y-<y>)


%% simulate classic full DCM (neural dynamic + measurement)
x = NaN;
while isweird(x)
    try
    x0;
    [y,x,x0,eta,e] = simulateNLSS(...
        n_t,f_fname,g_fname,theta,phi,u,alpha,sigma,options,x0);
    catch
    [y,x,x0,eta,e] = simulateNLSS(...
        n_t,f_fname,g_fname,theta,phi,u,alpha,sigma,options);    
    end
end

%% preparation for decoding simulation
try
    dim = options.dim;
catch
    dim.n_psi = length(psi);
    dim.n_t = n_t;
    dim.p = size(y,1);
end

dim.n = size(x0,1);
dim.n_r = options.inH.n_r;



%% initialize simul
% initial state

try 
    r0 ;
catch
    r0 = options.priors.muR0;
    sQ0 = getISqrtMat(options.priors.SigmaR0,0);
    r0 = r0 + sQ0*randn(dim.n_r,1);
end
% pre-allocate variables
r = zeros(dim.n_r,dim.n_t);


%% simulation
% Evaluate decoding function at initial conditions
if dim.n > 0
    %r(:,1) = VBA_evalFun('h',r0,psi,[x(:,1); u(:,1)],options,dim,1);

  
end

%[options,u,dim]= VBA_checkExtended(y,r,u,f_fname,g_fname,h_fname,dim,options);
% == move to VBA_Check
    options.h_fname=h_fname;
    options.h_nout=nargout(h_fname);
    options.checkGrads = 0;
% ==

% Evaluate decoding function at initial conditions
r(:,1) = VBA_evalFun('h',r0,psi,[x0 ; u(:,1)],options,dim,1);
for t = 2:dim.n_t
    r(:,t) = VBA_evalFun('h',r(:,t-1),psi,[x(:,t) ; u(:,t)],options,dim,t);
end
