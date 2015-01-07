function [y,x,x0,eta,e,u] = simulateNLSS(n_t,f_fname,g_fname,theta,phi,u,alpha,sigma,options,x0,fb)
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

% === Some verifications

try
    x0; 
    n=numel(x0);
catch
    try
        n = numel(options.priors.muX0);
    catch
        n = 0;
    end
end

feedback =  exist('fb','var');
if feedback
    try
    u(fb.indy,1);
    u(fb.indfb,1);
    catch
        error('*** Simulation with feedback requires allocated inputs');   
    end
end


% --- check dimensions
options = check_struct(options, ...
    'dim'   , struct, ...
    'priors', struct  ...
    ) ;
dim = options.dim;

dim = check_struct(dim, ...
    'n_theta'  , length(theta)  , ...
    'n_phi'    , length(phi)    , ...
    'n'        , n     , ...
    'n_t'      , n_t              ...
    );

try, options.inG; catch, options.inG = []; end
try, U = u(:,1);  catch, U = zeros(size(u,1),1); end
dim.p = size(feval(g_fname,zeros(dim.n,1),phi,U,options.inG),1);

options.dim=dim;
% --- check priors
options.priors = check_struct(options.priors, ...
    'a_alpha', 1, ...
    'b_alpha', 1  ...
    ) ;

if isinf(options.priors.a_alpha) && isequal(options.priors.b_alpha,0)
    options.priors.a_alpha = 1;
    options.priors.b_alpha = 1;
end

% --- check options
[options,u,dim] = VBA_check(zeros(dim.p,dim.n_t),u,f_fname,g_fname,dim,options);

% === Prepare simulation

% Get covariance structure
iQy = options.priors.iQy;
iQx = options.priors.iQx;

% Get time
et0 = clock;

% pre-allocate variables
x   = zeros(dim.n,dim.n_t);
eta = zeros(dim.n,dim.n_t);
e   = zeros(dim.p,dim.n_t);
y   = zeros(dim.p,dim.n_t);

% muxer
n_sources=numel(options.sources);
sgi = find([options.sources(:).type]==0) ;


% === Simulate timeseries

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
    x0 = zeros(dim.n,1);
end

% stack X0
x = [x0 x];

% Display progress
VBA_disp({ ...
    'Simulating SDE...'     , ...
	sprintf('%6.2f %%%%',0)     ...
 }, options);

%-- Loop over time points
for t = 1:dim.n_t
       
    % Evaluate evolution function at past hidden state   
    if dim.n > 0
        x(:,t+1) = VBA_evalFun('f',x(:,t),theta,u(:,t),options,dim,t) ;
        if ~isinf(alpha)
            Cx = VBA_getISqrtMat(iQx{t});
            eta(:,t) = (1./sqrt(alpha))*Cx*randn(dim.n,1);
            x(:,t+1) = x(:,t+1) + eta(:,t);
        end
    end

    % Evaluate observation function at current hidden state
    gt = VBA_evalFun('g',x(:,t+1),phi,u(:,t),options,dim,t);
      
    for i=1:n_sources
        s_idx = options.sources(i).out;
        switch options.sources(i).type
            % gaussian
            case 0 
                sigma_i = sigma(find(sgi==i)) ;
                if ~isinf(sigma_i)
                	C = VBA_getISqrtMat(iQy{t,find(sgi==i)});
                    e(s_idx,t) = (1./sqrt(sigma_i))*C*randn(length(s_idx),1);
                end
                y(s_idx,t) = gt(s_idx) + e(s_idx,t) ;
            % binomial
            case 1
                for k=1:length(s_idx)
                    y(s_idx(k),t) = sampleFromArbitraryP([gt(s_idx(k)),1-gt(s_idx(k))],[1,0]',1);
                end
        	% multinomial
            case 2
                resp = zeros(length(s_idx),1) ;
                resp(sampleFromArbitraryP(gt(s_idx),1:length(s_idx),1)) = 1;
                y(s_idx,t) = resp;
        end
        
    end
    e(:,t) = y(:,t) - gt;
    
    
    % fill in next input with last output and feedback
    if feedback && t < dim.n_t
        % get feedback on system's output
        if ~isempty(fb.h_fname)
            u(fb.indfb,t+1) = feval(fb.h_fname,y(:,t),t,fb.inH);
        end
        u(fb.indy,t+1) = y(:,t);
    end   
    
    % Display progress
    if mod(100*t/dim.n_t,10) <1 
        VBA_disp({ ...
            repmat('\b',1,9) ,  ...
            sprintf('%6.2f %%%%',floor(100*t/dim.n_t)), ...
        }, options);
    end
    
    if isweird({x(:,t),y(:,t)})
        break
    end
    
end

%unstack X0
x(:,1) = [];

% Display progress
VBA_disp({ ...
    repmat('\b',1,9)                                        ,...
	[' OK (took ',num2str(etime(clock,et0)),' seconds).']    ...
 },options);



