function [muy,Vy,iVp] = VBA_getLaplace(u,f_fname,g_fname,dim,options)

% returns the Laplace approximation to the prior predictive density
% function [muy,Vy,iVp] = VBA_getLaplace(u,f_fname,g_fname,dim,options)
% IN:
%   - u: experimentally controlled input (design)
%   - f_fname: the evolution function
%   - g_fname: the observation function
%   - dim: the model dimension structure
%   - options: the options structure
% OUT:
%   - muy: the 1st-order moment of the prior predictive density.
%   - Vy: the second-order moement of the prior predictive density. Note
%   that (time) lagged covariances are neglected. This gives a
%   block-diagonal structure to the data covariance matrix.
%   - iVp: the predicted posterior precision matrix of the model parameters

if nargout>2
    getiVp = 1;
else
    getiVp = 0;
end

options.checkGrads = 0; % well, this should have been done before...
options.priors.a_alpha = 0; % to bypass ODE transform in VBA_check.m
options.verbose = 0; % to quicken VBA_check.m
% try
%     delays = options.delays;
% catch
%     delays = [];
% end
% options.delays = [];
[options,u,dim] = VBA_check([],u,f_fname,g_fname,dim,options);

% Get prior covariance matrix
Sigma = zeros(dim.n_phi+dim.n_theta+dim.n,dim.n_phi+dim.n_theta+dim.n);
if dim.n_phi > 0
    options.priors.SigmaPhi
    Sigma(1:dim.n_phi,1:dim.n_phi)
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
if getiVp
    iVp = pinv(Sigma);
end

% Obtain derivatives of observations wrt...
dgdp = zeros(dim.n_phi+dim.n_theta+dim.n,dim.p);
x = zeros(dim.n,dim.n_t);
gx = zeros(dim.p,dim.n_t);

% initial condition
if dim.n > 0
    x0 = options.priors.muX0;
    Theta = options.priors.muTheta;
    [x(:,1),dF_dX,dF_dP] = ...
        VBA_evalFun('f',x0,Theta,u(:,1),options,dim,1);
    % get gradients wrt states
    dxdx0 = dF_dX;
    dxdTheta = dF_dP;
end
Phi = options.priors.muPhi;
[gx(:,1),dG_dX,dG_dP] = ...
    VBA_evalFun('g',x(:,1),Phi,u(:,1),options,dim,1);

% get gradients wrt to observations
if dim.n_phi > 0
    dgdp(1:dim.n_phi,:) = dG_dP;
end
if dim.n_theta > 0
    dgdp(dim.n_phi+1:dim.n_phi+dim.n_theta,:) = dxdTheta*dG_dX;
end
if dim.n > 0
    dgdp(dim.n_phi+dim.n_theta+1:end,:) = dxdx0*dG_dX;
end
muy(1:dim.p,1) = gx(:,1);
if options.binomial
    gx(:,1) = checkGX(gx(:,1)); % fix numerical instabilities
    Vy(1:dim.p,1:dim.p) = diag(gx(:,1).*(1-gx(:,1)));
    tmp = 1./gx(:,1) + 1./(1-gx(:,1));
    if getiVp
        iVp = iVp + dG_dP*diag(tmp)*dG_dP';
    end
else
    varY = options.priors.b_sigma./options.priors.a_sigma;
    Qy = VB_inv(options.priors.iQy{1});
    Vy(1:dim.p,1:dim.p) = dgdp'*Sigma*dgdp + varY.*Qy;
    if getiVp
        iVp = iVp + dgdp*options.priors.iQy{1}*dgdp'./varY;
    end
end

for t = 2:dim.n_t
    if dim.n > 0
        [x(:,t),dF_dX,dF_dP] = ...
            VBA_evalFun('f',x(:,t-1),Theta,u(:,t),options,dim,t);
        % Obtain derivatives of path wrt parameters...
         if dim.n_theta > 0
             dxdTheta = dF_dP + dxdTheta*dF_dX;
         end
        % ... and initial conditions
        dxdx0 = dxdx0*dF_dX;
    end
    [gx(:,t),dG_dX,dG_dP] = ...
        VBA_evalFun('g',x(:,t),Phi,u(:,t),options,dim,t);
    if dim.n_phi > 0
        dgdp(1:dim.n_phi,:) = dG_dP;
    end
    if dim.n_theta > 0
        dgdp(dim.n_phi+1:dim.n_phi+dim.n_theta,:) = dxdTheta*dG_dX;
    end
    if dim.n > 0
        dgdp(dim.n_phi+dim.n_theta+1:end,:) = dxdx0*dG_dX;
    end
    muy(dim.p*(t-1)+1:dim.p*t) = gx(:,t);
    if options.binomial
        gx(:,t) = checkGX(gx(:,t)); % fix numerical instabilities
        Vy(dim.p*(t-1)+1:dim.p*t,dim.p*(t-1)+1:dim.p*t) = ...
            diag(gx(:,t).*(1-gx(:,t)));
        if getiVp
            tmp = 1./gx(:,t) + 1./(1-gx(:,t));
            iVp = iVp + dG_dP*diag(tmp)*dG_dP';
        end
    else
        Qy = VB_inv(options.priors.iQy{t});
        Vy(dim.p*(t-1)+1:dim.p*t,dim.p*(t-1)+1:dim.p*t) = ...
            dgdp'*Sigma*dgdp + varY.*Qy;
        if getiVp
            iVp = iVp + dgdp*options.priors.iQy{t}*dgdp'./varY;
        end
    end
end

function x = checkGX(x)
x(x==0) = eps;
x(x==1) = 1-eps;


