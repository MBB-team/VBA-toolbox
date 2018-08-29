function [muy,Vy,iVp] = VBA_getLaplace(u,f_fname,g_fname,dim,options,checkVar,reduceVy)
% returns the Laplace approximation to the prior predictive density
% function [muy,Vy,iVp] = VBA_getLaplace(u,f_fname,g_fname,dim,options)
% IN:
%   - u: experimentally controlled input (design)
%   - f_fname: the evolution function
%   - g_fname: the observation function
%   - dim: the model dimension structure
%   - options: the options structure
%   - checkVar: flag for eyeballing the quality of the covariance matrix
% OUT:
%   - muy: the 1st-order moment of the prior predictive density.
%   - Vy: the second-order moment of the prior predictive density.
%   - iVp: the predicted posterior precision matrix of the model parameters
% SEE ALSO: splitLaplace

% _________________________________________________________________________
% checks
if ~exist('checkVar'), checkVar = 0; end
if ~exist('reduceVy'), reduceVy = 'full'; end

options.checkGrads     = 0; % well, this should have been done before...
options.priors.a_alpha = 0; % to bypass ODE transform in VBA_check.m
options.verbose        = 0; % to quicken VBA_check.m
[options,u,dim]        = VBA_check([],u,f_fname,g_fname,dim,options);

get_iVp = (nargout >= 3);

% _________________________________________________________________________
% initialization

% + memory preallocations
muy      = zeros(dim.p.*dim.n_t,1);                                     % first moment of the predictive density
switch reduceVy
    case 'full'
        Vy       = zeros(dim.p.*dim.n_t);                               % second moment of the predictive density
        dgdp     = zeros(dim.n_phi+dim.n_theta+dim.n,dim.p*dim.n_t);    % derivatives of obs wrt to paramters
    case 'diag'
        Vy       = zeros(dim.p.*dim.n_t,1);   
        dgdp     = zeros(dim.n_phi+dim.n_theta+dim.n,dim.p);            % derivatives of obs wrt to paramters
    case 'skip'
        Vy = NaN;
end

%x        = zeros(dim.n,1);                                             % current hidden state
%gx       = zeros(dim.p,1);                                             % current observation
dxdTheta = zeros(dim.n_theta,dim.n);                                    % derivatives of states wrt to theta
dxdX0    = eye(dim.n,dim.n);                                            % derivatives of states wrt to initial state
dF_dP    = nan(dim.n_theta,dim.n);                            
dF_dX    = nan(dim.n,dim.n);                            

% + initial values
Phi   = options.priors.muPhi;
Theta = options.priors.muTheta;
X0    = options.priors.muX0;

Sigma = blkdiag(options.priors.SigmaPhi         , ...
                options.priors.SigmaTheta       , ...
                options.priors.SigmaX0          );              % prior covariance matrix

if get_iVp
   iVp = VBA_inv(Sigma);                                        % predicted posterior precision matrix of the model parameters
end

gsi = find([options.sources(:).type]==0);

% _________________________________________________________________________
% integration loop

% start with initial state
x = X0;

% loop over time
for t = 1:dim.n_t
    
    % shortcuts
    offset = dim.p*(t-1) ;
    
    % .....................................................................
    % apply evolution and observation function
    if dim.n > 0
        [x,dF_dX,dF_dP] = VBA_evalFun('f',x,Theta,u(:,t),options,dim,t);  
    end
    [gx,dG_dX,dG_dP] = VBA_evalFun('g',x,Phi,u(:,t),options,dim,t);
    
    % .....................................................................
    % Obtain derivatives of path wrt parameters and initial conditions
    dxdTheta = dF_dP + dxdTheta*dF_dX;
    dxdX0 = dxdX0*dF_dX;

    % .....................................................................
    % store first moment and parameter gradient
    muy(offset+(1:dim.p)) = gx;
    dgdp_t = [dG_dP ; dxdTheta*dG_dX ; dxdX0*dG_dX];
    
    % .....................................................................
    % compute second moment
    for si=1:numel(options.sources)
        s_idx = options.sources(si).out ;
        s_idx = s_idx(options.isYout(s_idx,t)==0);
        
        if ~isempty(s_idx)
            
            s_idx_t = offset+s_idx ;
            
            switch options.sources(si).type
                case 0 % gaussian
                    g_idx = find(gsi==si);
                    varY = options.priors.b_sigma(g_idx)./options.priors.a_sigma(g_idx);
                    Qy = VBA_inv(options.priors.iQy{t,g_idx});
                    Vy_tmp = varY.*Qy;
                    tmp = options.priors.iQy{t,g_idx}./varY ;
                case 1 % binomial
                    gt = gx(s_idx);
                    Vy_tmp = diag(gt.*(1-gt));
                    tmp = diag(1./gt + 1./(1-gt));
                case 2 % multinomial
                    gt = gx(s_idx);
                    Vy_tmp = diag(gt.*(1-gt));
                    tmp = diag(1./gt);
            end
            
            % aggregate and store
            switch reduceVy
                case 'full'
                	Vy(s_idx_t,s_idx_t) = Vy_tmp ;    
                case 'diag'
                	Vy(s_idx_t) = Vy(s_idx_t) + diag(Vy_tmp) ;
                case 'skip'
            end

            if get_iVp
                iVp = iVp + dgdp_t(:,s_idx)*tmp*dgdp_t(:,s_idx)';
            end
        end
    end
    
    switch reduceVy
        case 'full'
            dgdp(:,offset+(1:dim.p)) = dgdp_t;
        case 'diag'
            Vy(offset+(1:dim.p)) = Vy(offset+(1:dim.p)) + diag(dgdp_t'*Sigma*dgdp_t);
        case 'skip'
    end

    
end

% _________________________________________________________________________
% form Laplace approximation to the covariance matrix
switch reduceVy
    case 'full'
        Vy = Vy + dgdp'*Sigma*dgdp;
    case 'diag'
    case 'skip'
end
            

% _________________________________________________________________________
% optional display

if checkVar
    in = struct('options',options);
    
    % get micro-time time series
    theta = options.priors.muTheta;
    phi = options.priors.muPhi;
    X0 = options.priors.muX0;
    
    % get prior covariance structure
    dgdtheta = VBA_numericDiff(@getObs,1,theta,phi,X0,u,in);
    dgdphi   = VBA_numericDiff(@getObs,2,theta,phi,X0,u,in);
    dgdX0    = VBA_numericDiff(@getObs,3,theta,phi,X0,u,in);
    
    Vy2 = dgdtheta'*options.priors.SigmaTheta*dgdtheta  ...
        + dgdphi'  *options.priors.SigmaPhi  *dgdphi    ...
        + dgdX0'   *options.priors.SigmaX0   *dgdX0     ...
        + Vy0;
    
    [hf] = VBA_displayGrads(Vy,Vy2);
end


function [gx] = getObs(theta,phi,X0,u,in)
in2 = struct('muTheta',theta,'muPhi',phi,'muX0',X0);
[x,gx,microTime,sampleInd] = VBA_microTime(in2,u,in);
gx = gx(:,sampleInd);
gx = gx(:);




