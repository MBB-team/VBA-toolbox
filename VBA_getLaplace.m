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

try; checkVar; catch; checkVar = 0; end

options.checkGrads = 0; % well, this should have been done before...
options.priors.a_alpha = 0; % to bypass ODE transform in VBA_check.m
options.verbose = 0; % to quicken VBA_check.m
[options,u,dim] = VBA_check([],u,f_fname,g_fname,dim,options);

% Get prior covariance matrix
Sigma = zeros(dim.n_phi+dim.n_theta+dim.n,dim.n_phi+dim.n_theta+dim.n);
if dim.n_phi > 0
    Sigma(1:dim.n_phi,1:dim.n_phi) = options.priors.SigmaPhi;
end
if dim.n_theta > 0
    Sigma(dim.n_phi+1:dim.n_phi+dim.n_theta,dim.n_phi+1:dim.n_phi+dim.n_theta) = options.priors.SigmaTheta;
end
if dim.n > 0
    Sigma(dim.n_phi+dim.n_theta+1:end,dim.n_phi+dim.n_theta+1:end) = options.priors.SigmaX0;
end

% pre-allocate output variables
Vy = zeros(dim.p.*dim.n_t);
muy = zeros(dim.p.*dim.n_t,1);
if nargout < 3
    get_iVp = 0;
else
    get_iVp = 1;
    iVp = VBA_inv(Sigma);
end


% Obtain derivatives of observations wrt...
dgdp = zeros(dim.n_phi+dim.n_theta+dim.n,dim.p*dim.n_t);
x = zeros(dim.n,dim.n_t);
gx = zeros(dim.p,dim.n_t);

% initial condition
if dim.n > 0
    x0 = options.priors.muX0;
    Theta = options.priors.muTheta;
    [x(:,1),dF_dX,dF_dP] = VBA_evalFun('f',x0,Theta,u(:,1),options,dim,1);
    % get gradients wrt states
    dxdx0 = dF_dX;
    dxdTheta = dF_dP;
end
Phi = options.priors.muPhi;
[gx(:,1),dG_dX,dG_dP] = VBA_evalFun('g',x(:,1),Phi,u(:,1),options,dim,1);

% get gradients wrt to observations
dgdp(1:dim.n_phi,1:dim.p) = dG_dP;
if dim.n_theta > 0
    dgdp(dim.n_phi+1:dim.n_phi+dim.n_theta,1:dim.p) = dxdTheta*dG_dX;
end
if dim.n > 0
    dgdp(dim.n_phi+dim.n_theta+1:end,1:dim.p) = dxdx0*dG_dX;
end
muy(1:dim.p,1) = gx(:,1);
gsi = find([options.sources(:).type]==0);
for si=1:numel(options.sources)
    s_idx = options.sources(si).out ;
    s_idx = s_idx(options.isYout(s_idx,1)==0);
    if ~isempty(s_idx)
    % - binary
    if options.sources(si).type 
        % true binomial
        if length(options.sources(si).out) == 1 
            Vy_tmp = diag(gx(s_idx,1).*(1-gx(s_idx,1)));
            tmp = diag(1./gx(s_idx,1) + 1./(1-gx(s_idx,1)));
        % multinomial
        else 
            Vy_tmp = diag(gx(s_idx,1).*(1-gx(s_idx,1)));
            tmp = diag(1./gx(s_idx,1)); % MUXER_TODO : check this !
        end
    % - gaussian
    else 
        varY = options.priors.b_sigma(gsi(si))./options.priors.a_sigma(gsi(si));
        Qy = VBA_inv(options.priors.iQy{1,si});
        Vy_tmp = varY.*Qy;
        tmp = options.priors.iQy{1,si}./varY ;
    end
    
    % aggregate
    Vy(s_idx,s_idx) = Vy_tmp ;
    if get_iVp
        iVp = iVp + dgdp(:,s_idx)*tmp*dgdp(:,s_idx)';
    end
    end
end

for t = 2:dim.n_t
    if dim.n > 0
        [x(:,t),dF_dX,dF_dP] = VBA_evalFun('f',x(:,t-1),Theta,u(:,t),options,dim,t);
        if dim.n_theta > 0
            % Obtain derivatives of path wrt parameters...
            dxdTheta = dF_dP + dxdTheta*dF_dX;
        end
        % ... and initial conditions
        dxdx0 = dxdx0*dF_dX;
    end
    [gx(:,t),dG_dX,dG_dP] = VBA_evalFun('g',x(:,t),Phi,u(:,t),options,dim,t);
    
    dgdp(1:dim.n_phi,1+(t-1)*dim.p:t*dim.p) = dG_dP;
    if dim.n_theta > 0
        dgdp(dim.n_phi+1:dim.n_phi+dim.n_theta,1+(t-1)*dim.p:t*dim.p) = dxdTheta*dG_dX;
    end
    if dim.n > 0
        dgdp(dim.n_phi+dim.n_theta+1:end,1+(t-1)*dim.p:t*dim.p) = dxdx0*dG_dX;
    end
    muy(dim.p*(t-1)+1:dim.p*t) = gx(:,t);

    for si=1:numel(options.sources)
        s_idx = options.sources(si).out ;
        s_idx = s_idx(options.isYout(s_idx,1)==0);
        if ~isempty(s_idx)
        s_idx_t = dim.p*(t-1)+s_idx ;
        % - binary
        if options.sources(si).type
            % true binomial
            if length(options.sources(si).out) == 1
                Vy_tmp = diag(gx(s_idx,t).*(1-gx(s_idx,t)));
                tmp = diag(1./gx(s_idx,t) + 1./(1-gx(s_idx,t)));
                % multinomial
            else
                Vy_tmp = diag(gx(s_idx,t).*(1-gx(s_idx,t)));
                tmp = diag(1./gx(s_idx,t)); % MUXER_TODO : check this !
            end
            % - gaussian
        else
            varY = options.priors.b_sigma(gsi(si))./options.priors.a_sigma(gsi(si));
            Qy = VBA_inv(options.priors.iQy{t,si});
            Vy_tmp = varY.*Qy;
            tmp = options.priors.iQy{t,si}./varY ;
        end
        
        % aggregate
        Vy(s_idx_t,s_idx_t) = Vy_tmp ;
        if get_iVp
            iVp = iVp + dgdp(:,s_idx_t)*tmp*dgdp(:,s_idx_t)';
        end
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
    % get prior covariance structure
    dgdtheta = numericDiff(@getObs,1,theta,phi,x0,u,in);
    dgdphi = numericDiff(@getObs,2,theta,phi,x0,u,in);
    dgdx0 = numericDiff(@getObs,3,theta,phi,x0,u,in);
    Vy2 = dgdtheta'*options.priors.SigmaTheta*dgdtheta ...
        + dgdphi'*options.priors.SigmaPhi*dgdphi ...
        + dgdx0'*options.priors.SigmaX0*dgdx0 ...
        + Vy0;
    [hf] = VBA_displayGrads(Vy,Vy2);
end


function [gx] = getObs(theta,phi,x0,u,in)
in2 = struct('muTheta',theta,'muPhi',phi,'muX0',x0);
[x,gx,microTime,sampleInd] = VBA_microTime(in2,u,in);
gx = gx(:,sampleInd);
gx = gx(:);



