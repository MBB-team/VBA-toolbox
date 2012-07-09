function [err] = VBA_check_errors_old(y,u,f_fname,g_fname,dim,options)
% at this point, all options and priors have been set to perform inversion.
% However, errors may still happen due to
% - bad model specification
% - bad choice of priors on model variables
% This script tests for such errors and returns a specific message to
% inform about the type of error.



%-- Errors of model specification


% Get time
et0 = clock;

% pre-allocate variables
x = zeros(dim.n,dim.n_t);
eta = zeros(dim.n,dim.n_t);
e = zeros(dim.p,dim.n_t);
y = zeros(dim.p,dim.n_t);

theta = options.priors.muTheta;
phi = options.priors.muPhi;

% Initial hidden-states value
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


% Evaluate evolution function at initial conditions
if dim.n > 0
    try
        x(:,1) = VBA_evalFun('f',x0,theta,u(:,1),options,dim,1);
    catch ME
    disp('ERROR : there is a problem in your evolution function');
        disp('- check model dimensions (hidden states and parameters)');
        disp('- check indices in the function');

    disp('----------------')      
    disp([ME.getReport]);
    disp('----------------') 
        
        err= 1;
        return;
    end
    C = getISqrtMat(iQx{1});
    eta(:,1) = (1./sqrt(alpha))*C*randn(dim.n,1);
    x(:,1) = x(:,1) + eta(:,1);
end

% Evaluate observation function at x(:,1)
try
    gt = VBA_evalFun('g',x(:,1),phi,u(:,1),options,dim,1);
catch ME
    disp('ERROR : there is a problem in your observation function');
            disp('- check model dimensions (hidden states and parameters)');
        disp('- check indices in the function');
    disp('----------------')      
    disp([ME.getReport]);
    disp('----------------')      

    err= 1;
    return;
    
end
if ~options.binomial
    C = getISqrtMat(iQy{1});
    e(:,1) = (1./sqrt(sigma))*C*randn(dim.p,1);
    y(:,1) = gt + e(:,1);
else
    for i=1:dim.p
        y(i,1) = sampleFromArbitraryP([gt(i),1-gt(i)],[1,0],1);
    end
    e(:,1) = y(:,1) - gt;
end

%-- Loop over time points

for t = 2:dim.n_t
    
    % Evaluate evolution function at past hidden state
    if dim.n > 0
        Cx = getISqrtMat(iQx{t});
        eta(:,t) = (1./sqrt(alpha))*Cx*randn(dim.n,1);
        try
            x(:,t) = VBA_evalFun('f',x(:,t-1),theta,u(:,t),options,dim,t) + eta(:,t);
        catch ME
    disp('ERROR : there is a problem in your evolution function');
            disp('- check model dimensions (hidden states and parameters)');
        disp('- check indices in the function');
    disp('----------------')      
    disp([ME.getReport]);
    disp('----------------')              
    err= 1;
            return;
        end
    end
    
    % Evaluate observation function at current hidden state
    try
        gt = VBA_evalFun('g',x(:,t),phi,u(:,t),options,dim,t);
    catch ME
        
            disp('ERROR : there is a problem in your observation function');
                    disp('- check model dimensions (hidden states and parameters)');
        disp('- check indices in the function');
    disp('----------------')      
    disp([ME.getReport]);
    disp('----------------') 
        err= 1;
        return;
        
    end
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
    
    % fill in next input with last output and feedback
    
    
    
    if isweird({x,y})
        break
    end
    err = 0;
end





