function [y,x,x0,eta,e,u] = VBA_simulate (n_t,f_fname,g_fname,theta,phi,u,alpha,sigma,options,x0,fb)
% samples times series from sDCM generative model
% [y,x,x0,dTime,eta,eta0] =
% VBA_simulate (n_t,f_fname,g_fname,theta,phi,u,alpha,sigma,options,x0)
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
%   - fb: an optional feedback struture that contains the following fields:
%       .h_fname: the name/handle of the function that implements the
%       feedback mapping, i.e. that maps the system's output y_t to its
%       feedback h(y_t,t,inH) 
%       .inH: an optional entry structure for the feedback mapping
%       .indy: the vector of indices that are used to address the previous
%       system output (y_t-1) within the current input (u_t)
%       .indfb: the vector of indices that are used to address the feedback
%       h(y_t-1,t-1,inH) to the previous system output within the current
%       input (u_t) 
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
options = VBA_check_struct(options, ...
    'dim'   , struct, ...
    'priors', struct  ...
    ) ;
dim = options.dim;

dim = VBA_check_struct(dim, ...
    'n_theta'  , length(theta)  , ...
    'n_phi'    , length(phi)    , ...
    'n'        , n     , ...
    'n_t'      , n_t              ...
    );

try, options.inG; catch, options.inG = struct (); end
try, U = u(:,1);  catch, U = zeros(size(u,1),1); end
dim.p = size(g_fname(zeros(dim.n,1),phi,U,options.inG),1);

options.dim = dim;
% --- check priors
options.priors = VBA_check_struct (options.priors, ...
    'a_alpha', 1, ...
    'b_alpha', 1  ...
    ) ;

if isempty (options.priors.a_alpha) || (isinf (options.priors.a_alpha) && isequal (options.priors.b_alpha, 0))
    options.priors.a_alpha = 1;
    options.priors.b_alpha = 1;
end

% --- check options
[options, u, dim] = VBA_check (zeros (dim.p, dim.n_t), u, f_fname, g_fname, dim, options);

% === Prepare simulation

% Get covariance structure
iQy = options.priors.iQy;
iQx = options.priors.iQx;

% Get time
et0 = clock;

% pre-allocate variables
x   = zeros (dim.n, dim.n_t);
eta = zeros (dim.n, dim.n_t);
e   = zeros (dim.p, dim.n_t);
y   = zeros (dim.p, dim.n_t);

% muxer
n_sources = numel (options.sources);
sgi = find ([options.sources(:).type] == 0) ;

% === Simulate timeseries

% Initial hidden-states value
if dim.n > 0
    try
        x0;
        assert(~ VBA_isWeird(x0));
    catch
        x0 = VBA_random ('Gaussian', options.priors.muX0, options.priors.SigmaX0);
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
 }, options, false);

%-- Loop over time points
for t = 1:dim.n_t
       
    % Evaluate evolution function at past hidden state   
    if dim.n > 0 
        % stochastic innovation
        Cx = VBA_inv (iQx{t}) / alpha ;
        eta(:, t) = VBA_random ('Gaussian', zeros (dim.n, 1), Cx) ;
        % evolution
        x(:, t + 1) = ... 
            VBA_evalFun ('f', x(:, t), theta, u(:, t), options, dim, t) ...
            + eta(:, t);
    end

    % Evaluate observation function at current hidden state
    gt = VBA_evalFun ('g', x(:, t + 1), phi, u(:, t), options, dim, t);
      
    for i = 1 : n_sources
        s_idx = options.sources(i).out;
        switch options.sources(i).type
            % gaussian
            case 0 
                C = VBA_inv (iQy{t, sgi == i}) / sigma(sgi == i) ;
                y(s_idx,t) =  VBA_random ('Gaussian', gt(s_idx), C) ;
            % binomial
            case 1
                y(s_idx, t) = VBA_random ('Bernoulli', gt(s_idx));
                
        	% multinomial
            case 2
                y(s_idx, t) = VBA_random ('Multinomial', 1, gt(s_idx));
        end
        
    end
    e(:,t) = y(:,t) - gt;
    
    % fill in next input with last output and feedback
    if feedback && t < dim.n_t
        % get feedback on system's output
        if ~ isempty (fb.h_fname)
            u(fb.indfb, t + 1) = fb.h_fname (y(:, t), t, fb.inH);
        end
        u(fb.indy, t + 1) = y(:, t);
    end   
    
    % Display progress
    if mod(100*t/dim.n_t,10) <1 
        VBA_disp({ ...
            repmat('\b',1,8) ,  ...
            sprintf('%6.2f %%%%',floor(100*t/dim.n_t)), ...
        }, options, false);
    end
    
end

%unstack X0
x(:, 1) = [];

% checks
if VBA_isWeird (x)
    error('VBA_simulate: evolution function produced weird values!');
end
        
% Display progress
VBA_disp({ ...
    repmat('\b',1,8)                                        ,...
	[' OK (took ',num2str(etime(clock,et0)),' seconds).']    ...
 },options);



