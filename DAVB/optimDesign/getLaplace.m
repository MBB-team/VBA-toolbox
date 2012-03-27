function [muy,Vy,iVp] = getLaplace(u,f_fname,g_fname,dim,options,checkVar)
% returns the Laplace approximation to the prior predictive density
% function [muy,Vy,iVp] = getLaplace(u,f_fname,g_fname,dim,options)
% IN:
%   - u: experimentally controlled input (design)
%   - f_fname: the evolution function
%   - g_fname: the observation function
%   - dim: the model dimension structure
%   - options: the options structure
%   - checkVar: flag for eyeballing the quaity of the covariance matrix
% OUT:
%   - muy: the 1st-order moment of the prior predictive density.
%   - Vy: the second-order moement of the prior predictive density.
%   - iVp: the predicted posterior precision matrix of the model parameters
% SEE ALSO: splitLaplace
%------------------------------------------------------------
% Copyright (C) 2012 Jean Daunizeau / License GNU GPL v2
%------------------------------------------------------------

try; checkVar; catch; checkVar = 0; end

% options.checkGrads = 0; % well, this should have been done before...
options.priors.a_alpha = 0; % to bypass ODE transform in VBA_check.m
options.verbose = 0; % to quicken VBA_check.m
[options,u,dim] = VBA_check([],u,f_fname,g_fname,dim,options);

% Get prior covariance matrix
Sigma = zeros(dim.n_phi+dim.n_theta+dim.n,dim.n_phi+dim.n_theta+dim.n);
if dim.n_phi > 0
    Sigma(1:dim.n_phi,1:dim.n_phi) = options.priors.SigmaPhi;
end
if dim.n_theta > 0
    Sigma(dim.n_phi+1:dim.n_phi+dim.n_theta,...
        dim.n_phi+1:dim.n_phi+dim.n_theta) = options.priors.SigmaTheta;
end
if dim.n > 0
    Sigma(dim.n_phi+dim.n_theta+1:end,dim.n_phi+dim.n_theta+1:end) = ...
        options.priors.SigmaX0;
end

% pre-allocate output variables
Vy = zeros(dim.p.*dim.n_t);
muy = zeros(dim.p.*dim.n_t,1);
if nargout < 3
    get_iVp = 0;
else
    get_iVp = 1;
    iVp = pinv(Sigma);
end


% Obtain derivatives of observations wrt...
dgdp = zeros(dim.n_phi+dim.n_theta+dim.n,dim.p*dim.n_t);
x = zeros(dim.n,dim.n_t);
gx = zeros(dim.p,dim.n_t);

% initial condition
if dim.n > 0
    x0 = options.priors.muX0;
    Theta = options.priors.muTheta;
    [x(:,1),dF_dX,dF_dP] = ...
        VBA_evalFun('f',x0,Theta,u(:,1),options,dim);
    % get gradients wrt states
    dxdx0 = dF_dX;
    dxdTheta = dF_dP;
end
Phi = options.priors.muPhi;
[gx(:,1),dG_dX,dG_dP] = ...
    VBA_evalFun('g',x(:,1),Phi,u(:,1),options,dim);

% get gradients wrt to observations
dgdp(1:dim.n_phi,1:dim.p) = dG_dP;
if dim.n_theta > 0
    dgdp(dim.n_phi+1:dim.n_phi+dim.n_theta,1:dim.p) = dxdTheta*dG_dX;
end
if dim.n > 0
    dgdp(dim.n_phi+dim.n_theta+1:end,1:dim.p) = dxdx0*dG_dX;
end
muy(1:dim.p,1) = gx(:,1);
if options.binomial
    gx(:,1) = checkGX(gx(:,1)); % fix numerical instabilities
    Vy(1:dim.p,1:dim.p) = diag(gx(:,1).*(1-gx(:,1)));
    tmp = 1./gx(:,1) + 1./(1-gx(:,1));
    if get_iVp
        iVp = iVp + dG_dP(:,1:dim.p)*diag(tmp)*dG_dP(:,1:dim.p)';
    end
else
    varY = options.priors.b_sigma./options.priors.a_sigma;
    Qy = pinv(options.priors.iQy{1});
    Vy(1:dim.p,1:dim.p) = varY.*Qy;
    if get_iVp
        iVp = iVp + ...
            dgdp(:,1:dim.p)*options.priors.iQy{1}*dgdp(:,1:dim.p)'./varY;
    end
end

for t = 2:dim.n_t
    if dim.n > 0
        [x(:,t),dF_dX,dF_dP] = ...
            VBA_evalFun('f',x(:,t-1),Theta,u(:,t),options,dim);
        if dim.n_theta > 0
            % Obtain derivatives of path wrt parameters...
            dxdTheta = dF_dP + dxdTheta*dF_dX;
        end
        % ... and initial conditions
        dxdx0 = dxdx0*dF_dX;
    end
    [gx(:,t),dG_dX,dG_dP] = ...
        VBA_evalFun('g',x(:,t),Phi,u(:,t),options,dim);
    
    dgdp(1:dim.n_phi,1+(t-1)*dim.p:t*dim.p) = dG_dP;
    if dim.n_theta > 0
        dgdp(dim.n_phi+1:dim.n_phi+dim.n_theta,1+(t-1)*dim.p:t*dim.p) = ...
            dxdTheta*dG_dX;
    end
    if dim.n > 0
        dgdp(dim.n_phi+dim.n_theta+1:end,1+(t-1)*dim.p:t*dim.p) = ...
            dxdx0*dG_dX;
    end
    muy(dim.p*(t-1)+1:dim.p*t) = gx(:,t);
    if options.binomial
        % fix numerical instabilities
        gx(:,t) = checkGX(gx(:,t));
        Vy(dim.p*(t-1)+1:dim.p*t,dim.p*(t-1)+1:dim.p*t) = ...
            diag(gx(:,t).*(1-gx(:,t)));
        tmp = 1./gx(:,t) + 1./(1-gx(:,t));
        if get_iVp
            iVp = iVp + ...
                dgdp(:,1+(t-1)*dim.p:t*dim.p)*diag(tmp)*dgdp(:,1+(t-1)*dim.p:t*dim.p)';
        end
    else
        Qy = pinv(options.priors.iQy{t});
        Vy(dim.p*(t-1)+1:dim.p*t,dim.p*(t-1)+1:dim.p*t) = varY.*Qy;
        if get_iVp
            iVp = iVp + ...
                dgdp(:,1+(t-1)*dim.p:t*dim.p)*options.priors.iQy{t}*dgdp(:,1+(t-1)*dim.p:t*dim.p)'./varY;
        end
    end
end
% form Laplace approximation to the covariance matrix
Vy0 = Vy;
Vy = Vy0 + dgdp'*Sigma*dgdp;

if checkVar
    
    in = struct('options',options);
    
    % get micro-time time series
    theta = options.priors.muTheta;
    phi = options.priors.muPhi;
    x0 = options.priors.muX0;
%     muy2 = getObs(theta,phi,x0,u,in);
    
    % get prior covariance structure
    dgdtheta = numericDiff(@getObs,1,theta,phi,x0,u,in);
    dgdphi = numericDiff(@getObs,2,theta,phi,x0,u,in);
    dgdx0 = numericDiff(@getObs,3,theta,phi,x0,u,in);
    
    Vy2 = dgdtheta'*options.priors.SigmaTheta*dgdtheta ...
        + dgdphi'*options.priors.SigmaPhi*dgdphi ...
        + dgdx0'*options.priors.SigmaX0*dgdx0 ...
        + Vy0;
    
    
    [hf] = VB_displayGrads(Vy,Vy2);
    pause
    
    
end


function [gx] = getObs(theta,phi,x0,u,in)

in2 = struct('muTheta',theta,...
    'muPhi',phi,...
    'muX0',x0);
[x,gx,microTime,sampleInd] = VBA_microTime(in2,u,in);
gx = gx(:,sampleInd);
gx = gx(:);

function x = checkGX(x)
x(x==0) = eps;
x(x==1) = 1-eps;


