function [y,x,x0,eta,e] = simulateNLSS(n_t,f_fname,g_fname,theta,phi,u,alpha,sigma,options,x0)
% samples times series from sDCM generative model
% [y,x,x0,dTime,eta,eta0] =
% simulateNLSS(n_t,f_fname,g_fname,theta,phi,u,alpha,sigma,options,x0)
%
% This function creates the time series of hidden-states and measurements
% under the following nonlinear state-space model:
%   x_t = f(x_t-1,Theta,u_t) + f_t
%   y_t = g(x_t,Phi,u_t) + e_t
% where f and g are the evolution and observation, respectively.
% IN:
%   - n_t: the number of time bins for the time series of hidden-states and
%   observations, i.e. the time indices satisfy: 1<= t < n_t
%   - f_fname/g_fname: evolution/observation function names.
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
%   - e: the pxt measurement errors (e:=y-<y>)

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
        dim.n_r = length(options.inF.r);
    catch
        dim.n_r=0;
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
[options,u,dim] = VBA_check(zeros(dim.p,dim.n_t),u,f_fname,g_fname,dim,options);

% Get covariance structure
iQy = options.priors.iQy;
iQx = options.priors.iQx;

% Get time
et0 = clock;

% pre-allocate variables
x = zeros(dim.n,dim.n_t);
eta = zeros(dim.n,dim.n_t);
e = zeros(dim.p,dim.n_t);
y = zeros(dim.p,dim.n_t);

% muxer
n_sources=numel(options.sources);

% Initial hidden-states value
if dim.n > 0
    try
        x0;
    catch
        x0 = options.priors.muX0;
        sQ0 = VBA_getISqrtMat(options.priors.SigmaX0,0);
        x0 = x0 + sQ0*randn(dim.n,1);
        clear sQ0
    end
else
    x0 = [];
end

% Evaluate evolution function at initial conditions
if dim.n > 0
    x(:,1) = VBA_evalFun('f',x0,theta,u(:,1),options,dim,1);
    if ~isinf(alpha)
        C = VBA_getISqrtMat(iQx{1});
        eta(:,1) = (1./sqrt(alpha))*C*randn(dim.n,1);
        x(:,1) = x(:,1) + eta(:,1);
    end
end

% Evaluate observation function at x(:,1)
gt = VBA_evalFun('g',x(:,1),phi,u(:,1),options,dim,1);

sgi = find([options.sources(:).type]==0) ;
for i=1:n_sources
    s_idx = options.sources(i).out;
    if options.sources(i).type == 0 % gaussian
        sigma_i = sigma(find(sgi==i)) ;
        if ~isinf(sigma_i)
            C = getISqrtMat(iQy{1,find(sgi==i)});
            e(s_idx,1) = (1./sqrt(sigma_i))*C*randn(length(s_idx),1);
        end
        y(s_idx,1) = gt(s_idx) ;
        y(s_idx,1) = y(s_idx,1) + e(s_idx,1);
    elseif  options.sources(i).type == 1% binary
     % true binomial
            for k=1:length(s_idx)
                y(s_idx(k),1) = sampleFromArbitraryP([gt(s_idx(k)),1-gt(s_idx(k))],[1,0]',1);
                e(s_idx(k),1) = y(s_idx(k),1) - gt(s_idx(k));
            end
        else % multinomial
            resp = zeros(length(s_idx),1) ;
            resp(sampleFromArbitraryP(gt(s_idx),1:length(s_idx),1)) = 1;
            y(s_idx,1) = resp;
            e(s_idx,1) = y(s_idx,1) - gt(s_idx);
        
    end
    e(:,1) = y(:,1) - gt;
end

%-- Loop over time points

% Display progress
if options.verbose
    fprintf(1,'Simulating SDE...')
    fprintf(1,'%6.2f %%',0)
end

for t = 2:dim.n_t
        
    % Evaluate evolution function at past hidden state
    if dim.n > 0
        if ~isinf(alpha)
            Cx = VBA_getISqrtMat(iQx{t});
            eta(:,t) = (1./sqrt(alpha))*Cx*randn(dim.n,1);
        end
        x(:,t) = VBA_evalFun('f',x(:,t-1),theta,u(:,t),options,dim,t) + eta(:,t);
    end
    
    % Evaluate observation function at current hidden state
    gt = VBA_evalFun('g',x(:,t),phi,u(:,t),options,dim,t);
    
  
for i=1:n_sources
    s_idx = options.sources(i).out;
    if options.sources(i).type == 0 % gaussian
         sigma_i = sigma(find(sgi==i)) ;
        if ~isinf(sigma_i)
            C = VBA_getISqrtMat(iQy{t,find(sgi==i)});
            e(s_idx,t) = (1./sqrt(sigma_i))*C*randn(length(s_idx),1);
        end
        y(s_idx,t) = gt(s_idx) ;
        y(s_idx,t) = y(s_idx,t) + e(s_idx,t);
    else % binary
        if length(s_idx) == 1 % true binomial
%             y(s_idx,t) = sampleFromArbitraryP([gt(s_idx),1-gt(s_idx)],[1,0],1);
            y(s_idx,t) = binomial_sample(gt(s_idx));
        else % multinomial
            resp = zeros(length(s_idx),1) ;
            try
            resp(sampleFromArbitraryP(gt(s_idx),1:length(s_idx),1)) = 1;
            catch
                end
            y(s_idx,t) = resp;
        end
        e(s_idx,t) = y(s_idx,t) - gt(s_idx);
    end
end
    
    % Display progress
    if mod(100*t/dim.n_t,10) <1 && options.verbose
        fprintf(1,repmat('\b',1,8))
        fprintf(1,'%6.2f %%',floor(100*t/dim.n_t))
    end
    if isweird({x(:,t),y(:,t)})
        break
    end
    
end

% Display progress
if options.verbose
    fprintf(1,repmat('\b',1,8))
    fprintf(1,[' OK (took ',num2str(etime(clock,et0)),' seconds).'])
    fprintf(1,'\n')
end



