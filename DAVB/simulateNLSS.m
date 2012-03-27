function [y,x,x0,eta,e] = simulateNLSS(n_t,f_fname,g_fname,theta,phi,u,alpha,sigma,options,x0)
% samples times series from sDCM generative model
% [y,x,x0,dTime,eta,eta0] =
% simulateNLSS(n_t,f_fname,g_fname,theta,phi,u,alpha,sigma,options,x0)
%
% This function creates the time series of hidden-states and measurements
% under the following nonlinear state-space model:
%   y_t   = g(x_t,Phi,u_t,t) + e_t
%   x_t+1 = f(x_t,Theta,u_t,t) + f_t
% where f and g are the evolution and observation, respectively.
% IN:
%   - n_t: the number of time bins for the time series of hidden-states and
%   observations, i.e. the time indices satisfy: 1<= t < n_t
%   - f_fname/g_fname: evolution/observation function names. The time entry
%   of these functions IS NOT the time index, it is defined as :
%   t0 + t*delta_t.
%   - theta/phi: evolution/observation parameters values.
%   - u: the mxt input matrix
%   - alpha: precision of the stochastic innovations
%   - sigma: precision of the measurement error
%   - options: structure variable containing the following fields:
%       .inF
%       .inG
%   - x0: the initial conditions
% OUT:
%   - y: the pxt (noisy) measurement time series
%   - x: the nxt (noisy) hidden-states time series
%   - x0: the nx1 initial conditions
%   - eta: the nxt stochastic innovations time series
%   - e: the pxt measurement errors (for non binomial data) or the
%   likelihood p(y=1|P,m) for binomial data.
%------------------------------------------------------------
% Copyright (C) 2012 Jean Daunizeau / License GNU GPL v2
%------------------------------------------------------------
% 27/03/2007: JD.

% get system dimensions
try
    dim = options.dim;
catch
    dim.n_theta = length(theta);
    dim.n_phi = length(phi);
    dim.n_t = n_t;
    try
        dim.n = size(x0,1);
    catch
        dim.n = size(options.priors.muX0,1);
    end
    try
        options.inG;
    catch
        options.inG = [];
    end
    try
        U = u(:,1);
    catch
        U = [];
    end
    dim.p = size(feval(g_fname,zeros(dim.n,1),phi,U,options.inG),1);
end

try
    if isinf(options.priors.a_alpha) && isequal(options.priors.b_alpha,0)
        options.priors.a_alpha = 1;
        options.priors.b_alpha = 1;
    end
end
options.priors.AR = 0;
[options,u,dim] = VBA_check(zeros(dim.p,dim.n_t),u,f_fname,g_fname,dim,options);

% Get covariance structure
iQy = options.priors.iQy;
iQx = options.priors.iQx;

% Get time
et0 = clock;

%-- initial hidden-states value
if dim.n > 0
    try
        x0;
    catch
        x0 = options.priors.muX0;
        sQ0 = getISqrtMat(options.priors.SigmaX0,0);
        x0 = x0 + sQ0*randn(dim.n,1);
        clear sQ0
    end
else
    x0 = [];
end

% pre-allocate variables
x = zeros(dim.n,dim.n_t);
eta = zeros(dim.n,dim.n_t);
e = zeros(dim.p,dim.n_t);
y = zeros(dim.p,dim.n_t);
if dim.n > 0
    % Evaluate evolution function at initial conditions
    x(:,1) = VBA_evalFun('f',x0,theta,u(:,1),options,dim);
    C = getISqrtMat(iQx{1});
    eta(:,1) = (1./sqrt(alpha))*C*randn(dim.n,1);
    x(:,1) = x(:,1) + eta(:,1);
end
% Evaluates observation function at x(:,1)
gt = VBA_evalFun('g',x(:,1),phi,u(:,1),options,dim);
if ~options.binomial
    C = getISqrtMat(iQy{1});
    e(:,1) = (1./sqrt(sigma))*C*randn(dim.p,1);
    y(:,1) = gt + e(:,1);
else
    for i=1:dim.p
        y(i,1) = sampleFromArbitraryP([gt(i),1-gt(i)],[1,0],1);
    end
    e(:,1) = gt;
end

%-- Loop over time points
% Display progress
if options.verbose
    fprintf(1,'Simulating SDE...')
    fprintf(1,'%6.2f %%',0)
end
for t = 2:dim.n_t
    if dim.n > 0
        % Evaluate evolution function at past hidden state
        Cx = getISqrtMat(iQx{t});
        eta(:,t) = (1./sqrt(alpha))*Cx*randn(dim.n,1);
        x(:,t) = VBA_evalFun('f',x(:,t-1),theta,u(:,t),options,dim) + eta(:,t);
    end
    % Evaluate observation function at current hidden state
    gt = VBA_evalFun('g',x(:,t),phi,u(:,t),options,dim);
    if ~options.binomial
        Cy = getISqrtMat(iQy{t});
        e(:,t) = (1./sqrt(sigma))*Cy*randn(dim.p,1);
        y(:,t) = gt + e(:,t);
    else
        for i=1:dim.p
            y(i,t) = sampleFromArbitraryP([gt(i),1-gt(i)],[1,0],1);
        end
        e(:,t) = y(:,t) - gt;
    end
    % Display progress
    if mod(100*t/dim.n_t,10) <1 && options.verbose
        fprintf(1,repmat('\b',1,8))
        fprintf(1,'%6.2f %%',floor(100*t/dim.n_t))
    end
    if isweird({x,y})
        break
    end
end
% Display progress
if options.verbose
    fprintf(1,repmat('\b',1,8))
    fprintf(1,[' OK (took ',num2str(etime(clock,et0)),' seconds).'])
    fprintf(1,'\n')
end

