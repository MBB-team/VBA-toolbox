function [muy,Vy,iVp] = VBA_getLaplace(u,f_fname,g_fname,dim,options,checkVar)
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

options.checkGrads     = 0; % well, this should have been done before...
options.priors.a_alpha = 0; % to bypass ODE transform in VBA_check.m
options.verbose        = 0; % to quicken VBA_check.m
[options,u,dim]        = VBA_check([],u,f_fname,g_fname,dim,options);

get_iVp = (nargout >= 3);

% _________________________________________________________________________
% initialization

% + memory preallocations
muy      = zeros(dim.p.*dim.n_t,1);                             % first moment of the predictive density
Vy       = zeros(dim.p.*dim.n_t);                               % first moment of the predictive density
x        = zeros(dim.n,dim.n_t);                                % hidden state time-series
gx       = zeros(dim.p,dim.n_t);                                % observation time-series
dxdTheta = zeros(dim.n_theta,dim.n);                            % derivatives of states wrt to theta
dxdX0    = eye(dim.n,dim.n);                                    % derivatives of states wrt to initial state
dF_dP    = nan(dim.n_theta,dim.n);                            
dF_dX    = nan(dim.n,dim.n);                            
dgdp     = zeros(dim.n_phi+dim.n_theta+dim.n,dim.p*dim.n_t);    % derivatives of obs wrt to paramters

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

% append initial state
x = [X0 x];

% loop over time
for t = 1:dim.n_t
    
    % shortcuts
    offset = dim.p*(t-1) ;
    
    % .....................................................................
    % apply evolution and observation function
    if dim.n > 0
        [x(:,t+1),dF_dX,dF_dP] = VBA_evalFun('f',x(:,t),Theta,u(:,t),options,dim,t);  
    end
    [gx(:,t),dG_dX,dG_dP] = VBA_evalFun('g',x(:,t+1),Phi,u(:,t),options,dim,t);
    
    % .....................................................................
    % Obtain derivatives of path wrt parameters and initial conditions
    dxdTheta = dF_dP + dxdTheta*dF_dX;
    dxdX0 = dxdX0*dF_dX;

    % .....................................................................
    % store first moment and parameter gradient
    muy(offset+(1:dim.p)) = gx(:,t);
    dgdp(:,offset+(1:dim.p)) = [dG_dP ; dxdTheta*dG_dX ; dxdX0*dG_dX];
    
    % .....................................................................
    % compute second moment
    for si=1:numel(options.sources)
        s_idx = options.sources(si).out ;
        s_idx = s_idx(options.isYout(s_idx,t)==0);
        
        if ~isempty(s_idx)
            
            s_idx_t = offset+s_idx ;
            
            switch options.sources(si).type
                case 0 % gaussian
                    varY = options.priors.b_sigma(gsi==si)./options.priors.a_sigma(gsi==si);
                    Qy = VBA_inv(options.priors.iQy{t,si});
                    Vy_tmp = varY.*Qy;
                    tmp = options.priors.iQy{t,si}./varY ;
                case 1 % binomial
                    gt = gx(s_idx,t);
                    Vy_tmp = diag(gt.*(1-gt));
                    tmp = diag(1./gt + 1./(1-gt));
                case 2 % multinomial
                    gt = gx(s_idx,t);
                    Vy_tmp = diag(gt.*(1-gt));
                    tmp = diag(1./gt);
            end
            
            % aggregate and store
            Vy(s_idx_t,s_idx_t) = Vy_tmp ;
            if get_iVp
                iVp = iVp + dgdp(:,s_idx_t)*tmp*dgdp(:,s_idx_t)';
            end
        end
    end
    
end

% _________________________________________________________________________
% form Laplace approximation to the covariance matrix
Vy = Vy + dgdp'*Sigma*dgdp;

% _________________________________________________________________________
% optional display

if checkVar
    in = struct('options',options);
    
    % get micro-time time series
    theta = options.priors.muTheta;
    phi = options.priors.muPhi;
    X0 = options.priors.muX0;
    
    % get prior covariance structure
    dgdtheta = numericDiff(@getObs,1,theta,phi,X0,u,in);
    dgdphi   = numericDiff(@getObs,2,theta,phi,X0,u,in);
    dgdX0    = numericDiff(@getObs,3,theta,phi,X0,u,in);
    
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




