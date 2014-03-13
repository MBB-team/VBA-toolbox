function [IX,SigmaX,deltaMuX,suffStat] = VBA_IX(...
    X,y,posterior,suffStat,dim,u,Te,options)
% Laplace-Kalman smoother (Gauss-Newton update of hidden states)

%  Look-up which hidden states to update
indIn = options.params2update.x;

% Get precision parameters
alphaHat = posterior.a_alpha./posterior.b_alpha;
sigmaHat = posterior.a_sigma./posterior.b_sigma;
iQx = options.priors.iQx;
iQy = options.priors.iQy;

% Preallocate intermediate variables
In = speye(dim.n);
mStar = zeros(dim.n,dim.n_t);
m = zeros(dim.n,dim.n_t);
ohmega = zeros(dim.n,dim.n_t);
muTilde = zeros(dim.n,dim.n_t);
muX = zeros(dim.n,dim.n_t);
fx = zeros(dim.n,dim.n_t);
dF_dX = cell(dim.n_t,1);
gx = zeros(dim.p,dim.n_t);
dG_dX = cell(dim.n_t,1);
iB = cell(dim.n_t,1);
iR = cell(dim.n_t,1);
C = cell(dim.n_t,1);
SigmaX.current = cell(dim.n_t,1);
SigmaX.inter = cell(dim.n_t-1,1);
dy = zeros(dim.p,dim.n_t);
vy = zeros(dim.p,dim.n_t);
gx = zeros(dim.p,dim.n_t);
dy2 = 0;
dx = zeros(dim.n,dim.n_t);
fx = zeros(dim.n,dim.n_t-1);
dx2 = 0;
Sphid2gdphi2 = 0;
Sthetad2fdtheta2 = 0;
SXd2fdx2 = 0;
SXtdfdx = 0;
SXd2gdx2 = 0;
trSx = 0;
div = 0;

%--------------------- Forward pass -----------------------%
set(options.display.hm(1),'string',...
    'VB Gauss-Newton EKF: forward pass... ');
set(options.display.hm(2),'string','0%');
drawnow

%---- Initial condition ----%
% !! Attn: MF approx. on initial condition !!

% evaluate evolution function at current mode
if ~options.priors.AR
    [fx0,dF_dX0,dF_dTheta0] =...
        VBA_evalFun('f',posterior.muX0,posterior.muTheta,u(:,1),options,dim);
else
    clear VBA_smoothNLSS
    [fx0,dF_dX0,dF_dTheta0] =...
        VBA_evalFun('f',zeros(dim.n,1),posterior.muTheta,u(:,1),options,dim);
end

% evaluate observation function at current mode
[gx(:,1),dG_dX{1},dG_dPhi,d2G_dXdPhi] =...
    VBA_evalFun('g',X(:,1),posterior.muPhi,u(:,1),options,dim);

% check infinite precision transition pdf
iQ = VB_inv(iQx{1},indIn{1},'replace');
Q = VB_inv(iQ,indIn{1});
IN = diag(~~diag(iQ));

% mean-field terms
if ~options.ignoreMF
    Sphid2gdphi2 = Sphid2gdphi2 ...
        + trace(dG_dPhi*iQy{1}*dG_dPhi'*posterior.SigmaPhi);
    Sthetad2fdtheta2 = Sthetad2fdtheta2 ...
        + trace(dF_dTheta0*iQ*dF_dTheta0'*posterior.SigmaTheta);
end
SXd2gdx2 = SXd2gdx2 ...
    +trace(dG_dX{1}*iQy{1}*dG_dX{1}'*posterior.SigmaX.current{1});
SXd2fdx2 = SXd2fdx2 ...
    +trace(dF_dX0*iQ*dF_dX0'*posterior.SigmaX0);
trSx = trSx + trace(posterior.SigmaX.current{1});

% error terms
dx(:,1) = IN*(X(:,1) - fx0);
dx2 = dx2 + dx(:,1)'*iQ*dx(:,1);
dy(:,1) = y(:,1) - gx(:,1);
dy2 = dy2 + dy(:,1)'*iQy{1}*dy(:,1);

% covariance matrices
Rp = dF_dX0'*posterior.SigmaX0*dF_dX0 + (Te./alphaHat)*Q;
C{1} = (sigmaHat./Te)*dG_dX{1}*iQy{1}*dG_dX{1}';
if ~options.ignoreMF && dim.n_phi > 0
    A1g = reshape(permute(d2G_dXdPhi,[2,3,1]),dim.p*dim.n_phi,dim.n)';
    A2g = A1g*kron(iQy{1},posterior.SigmaPhi);
    C{1} = C{1} + (sigmaHat./Te)*A2g*A1g';
end
iRp{1} = VB_inv(Rp,indIn{1});
iR{1} =  iRp{1} + C{1};
R = VB_inv(iR{1},indIn{1});

% prediction : p(X_1|X_0)
mStar(:,1) = fx0;

% update : p(X_1|Y_1)
m(:,1) = mStar(:,1) ...
    + (sigmaHat./Te).*IN*R*dG_dX{1}*iQy{1}*...
    ( dy(:,1) + dG_dX{1}'*(X(:,1)-mStar(:,1)) );
if ~options.ignoreMF && dim.n_phi > 0
    m(:,1) = m(:,1) ...
        + (sigmaHat./Te).*IN*R*A2g*( A1g'*(X(:,1)-mStar(:,1)) -vec(dG_dPhi) );
end

% Transient posterior
SigmaX.current{1} = R;
muX(:,1) = m(:,1);

% Predictive density (data space)
V = (Te./sigmaHat).*pinv(iQy{1}) + dG_dX{1}'*R*dG_dX{1};
if dim.n_phi > 0
    V = V + dG_dPhi'*posterior.SigmaPhi*dG_dPhi;
end
vy(:,1) = diag(V);

% Entropy calculus
SX = 0.5*dim.n*dim.n_t*log(2*pi*exp(1)) ...
    + 0.5*VBA_logDet(posterior.SigmaX.current{1},indIn{1});

% % artificial prior
% muxp = zeros(dim.n,dim.n_t);
% SUxp = cell(dim.n_t,1);
% CUxp = cell(dim.n_t,1);
% muxp(:,1) = options.priors.muX(:,1);
% SUxp{1} = options.priors.SigmaX.current{1};


%---- Sequential message-passing algorithm: forward pass (filter) ----%
for t = 2:dim.n_t
   
    % check infinite precision transition pdf
    iQ = VB_inv(iQx{t},indIn{t},'replace');
    Q =  VB_inv(iQx{t},indIn{t});
    IN = diag(~~diag(iQ));
    
    % evaluate evolution function at current mode
    [fx(:,t-1),dF_dX{t-1},dF_dTheta,d2F_dXdTheta] =...
        VBA_evalFun('f',X(:,t-1),posterior.muTheta,u(:,t),options,dim);

    % evaluate observation function at current mode
    [gx(:,t),dG_dX{t},dG_dPhi,d2G_dXdPhi] =...
        VBA_evalFun('g',X(:,t),posterior.muPhi,u(:,t),options,dim);

%     % artificial prior
%     args = {'f',muxp(:,t-1),posterior.muTheta,u(:,t),options,dim};
%     [muxp(:,t),SUxp{t},CUxp{t}] = VB_UT(muxp(:,t-1),SUxp{t-1},'VBA_evalFun',args,2);
%     SUxp{t} = SUxp{t} + Q./alphaHat;
    
    
    % mean-field terms
    try
        Sphid2gdphi2 = Sphid2gdphi2 ...
            + trace(dG_dPhi*iQy{t}*dG_dPhi'*posterior.SigmaPhi);
    end
    try
        Sthetad2fdtheta2 = Sthetad2fdtheta2 ...
            + trace(dF_dTheta*iQ*dF_dTheta'*posterior.SigmaTheta);
    end
    SXd2gdx2 = SXd2gdx2 ...
        + trace(dG_dX{t}*iQy{t}*dG_dX{t}'*posterior.SigmaX.current{t});
    SXd2fdx2 = SXd2fdx2 ...
        + trace(dF_dX{t-1}*iQ*dF_dX{t-1}'*posterior.SigmaX.current{t});
    SXtdfdx = SXtdfdx ...
        - 2*trace( iQ*dF_dX{t-1}'*posterior.SigmaX.inter{t-1} );
    trSx = trSx + trace(posterior.SigmaX.current{t});

    % error terms
    dx(:,t) = (X(:,t) - fx(:,t-1));
    dx2 = dx2 + dx(:,t)'*iQ*dx(:,t);
    dy(:,t) = y(:,t) - gx(:,t);
    dy2 = dy2 + dy(:,t)'*iQy{t}*dy(:,t);

    % covariance matrices
    B = iR{t-1} + (alphaHat./Te)*dF_dX{t-1}*iQ*dF_dX{t-1}';
    if ~options.ignoreMF && dim.n_theta > 0
        A1f = reshape(permute...
            (d2F_dXdTheta,[2,3,1]),dim.n*dim.n_theta,dim.n)';
        A2f = A1f*kron(iQ,posterior.SigmaTheta);
        B = B + (alphaHat./Te)*A2f*A1f';
    end
    iB{t-1} = VB_inv(B,indIn{t});
    iRp2 = (alphaHat./Te)*iQ*...
        ( In - (alphaHat./Te)*iQ*dF_dX{t-1}'*iB{t-1}*dF_dX{t-1} );
    C{t} = (sigmaHat./Te)*dG_dX{t}*iQy{t}*dG_dX{t}';
    if ~options.ignoreMF && dim.n_phi > 0
        A1g = reshape(permute(d2G_dXdPhi,[2,3,1]),dim.p*dim.n_phi,dim.n)';
        A2g = A1g*kron(iQy{t},posterior.SigmaPhi);
        C{t} = C{t} + (sigmaHat./Te)*A2g*A1g';
    end
    
    Rp = dF_dX{t-1}'*R*dF_dX{t-1} + (Te./alphaHat)*Q;
    iRp{t} = VB_inv(Rp,indIn{t});
    iR{t} = iRp{t} + C{t};
    R = VB_inv(iR{t},indIn{t});
    
    % prediction : p(X_t|Y_1:t-1)
    mStar(:,t) = fx(:,t-1) ...
        + IN*dF_dX{t-1}'*(m(:,t-1) -X(:,t-1));
    if ~options.ignoreMF && dim.n_theta > 0
        mStar(:,t) = mStar(:,t) ...
            + IN*dF_dX{t-1}'*pinv((Te/alphaHat)*iR{t-1}+A2f*A1f')*A2f*...
            ( A1f'*(X(:,t-1)-m(:,t-1)) -vec(dF_dTheta) );
    end       
    
    % update : p(X_t|Y_1:t)
    m(:,t) = mStar(:,t) ...
        + (sigmaHat./Te).*IN*R*dG_dX{t}*iQy{1}*...
        ( dy(:,t) + dG_dX{t}'*(X(:,t)-mStar(:,t)) );
    if ~options.ignoreMF && dim.n_phi > 0
        m(:,t) = m(:,t) ...
            + (sigmaHat./Te).*IN*R*A2g*...
            ( A1g'*(X(:,t)-mStar(:,t)) -vec(dG_dPhi) );
    end
    
    % Transient posterior
    SigmaX.current{t} = R;
    muX(:,t) = m(:,t);

    % Predictive density (data space)
    V = (Te./sigmaHat).*pinv(iQy{t}) + dG_dX{t}'*R*dG_dX{t};
    if dim.n_phi > 0
        V = V + dG_dPhi'*posterior.SigmaPhi*dG_dPhi;
    end
    vy(:,t) = diag(V);

    % Display progress
    if mod(t,dim.n_t./10) == 0
        set(options.display.hm(2),'string',[num2str(100*t/dim.n_t),'%']);
        drawnow
    end
    
    % Accelerate divergent update
    if any(isinf(dy2)) || any(isnan(dy2)) || ...
            any(isinf(dx2)) || any(isnan(dx2)) || ...
            any(isinf(dG_dX{t}(:))) || any(isnan(dG_dX{t}(:))) || ...
            any(isinf(dF_dX{t-1}(:))) || any(isnan(dF_dX{t-1}(:))) || ...
            any(~isreal(dG_dX{t}(:))) || any(~isreal(dF_dX{t-1}(:)))
        div = 1;
        break
    end
    
end         %--- end of forward pass ---%

% Display progress
set(options.display.hm(2),'string','OK');
drawnow


%------------------------ backward pass ------------------------%
set(options.display.hm(1),'string',...
    'VB Gauss-Newton EKF: backward pass... ');
set(options.display.hm(2),'string','0%');
drawnow


%---- Sequential message-passing algorithm: backward pass (smoother) ----%
if ~div % only if forward-pass converged
    
    % check infinite precision transition pdf
    iQ = VB_inv(iQx{dim.n_t},indIn{dim.n_t},'replace');
    
    % Final condition beta-message (by convention)
%     ohmega(:,dim.n_t) = m(:,dim.n_t);
%     iS = iR{dim.n_t};
    E = (alphaHat./Te)*iQ + (sigmaHat./Te)*C{dim.n_t};

    for i = 1:dim.n_t-1

        % time index variable change
        t = dim.n_t-i;
        
        % check infinite precision transition pdf
        iQ = VB_inv(iQx{t},indIn{t},'replace');
        Qi = VB_inv(iQx{t},indIn{t});
        Qy = VB_inv(iQy{t+1},[]);
        iQ2 = VB_inv(iQx{t+1},indIn{t+1},'replace');
        IN = diag(~~diag(iQ));
        
%         % backward Markov kernel
%         tmp = CUxp{t+1}*VB_inv(SUxp{t+1},[]);
%         mub = muxp(:,t) + tmp*(X(:,t+1)-muxp(:,t+1));
%         Pb = SUxp{t} - tmp*CUxp{t+1}';
%         
%         % observational moments
%         args = {'g',X(:,t),posterior.muPhi,u(:,t),options,dim};
%         [muy,Py,Pxy] = VB_UT(mub,Pb,'VBA_evalFun',args,2);
%         Py = Py + Qy./sigmaHat;
%         
%         % backward update step
%         tmp = Pxy*VB_inv(Py,[]);
%         ohmega = mub + tmp*(y(:,t)-muy);
%         Px = Pb - tmp*Pxy';
%         
%         iP = iRp{t} + VB_inv(Px,[]);
%         xTilde = VB_inv(iP,[])*(iRp{t}*m(:,t) + iP*ohmega);
%         
%         SigmaX.current{t} = VB_inv(iP-VB_inv(SUxp{t},[]),[]);
%         muX(:,t) = SigmaX.current{t}*(iP*xTilde-SUxp{t}*muxp(:,t));
        
%         % call unscented transform
%         mX = [X(:,t);zeros(dim.n,1);zeros(dim.p,1)];
%         PX = sparse(2*dim.n:dim.p,2*dim.n:dim.p);
%         PX(1:dim.n,1:dim.n) = posterior.Sigma.current{t};
%         PX(dim.n+1:2*dim.n,dim.n+1:2*dim.n) = Qi./alphaHat;
%         PX(2*dim.n+1:2*dim.n+dim.p,2*dim.n+1:2*dim.n+dim.p) = Qy./sigmaHat;
%         args = {mX,posterior,u(:,t+1),options,dim};
%         [muU,SU,CU] = VB_UT(mX,PX,'getGFX',args,1);
%         
%         mux = muU(1:dim.n);
%         muy = muU(2*dim.n+1:2*dim.n+dim.p);
%         Py = SU;
%         Pxy = 
        
%         % call unscented transform
%         args = {'f',X(:,t),posterior.muTheta,u(:,t+1),options,dim};
%         [muU,SU,CU] = VB_UT(X(:,t),posterior.SigmaX.current{t},'VBA_evalFun',args,2);
%         SU = SU + Qi./alphaHat;
%         
%         D = CU*VB_inv(SU,[]);
%         SigmaX.current{t} = SigmaX.current{t} + D*(SigmaX.current{t+1} - SU)*D';
%         muX(:,t) = m(:,t) + D*(muX(:,t+1)-muU);
        
        
        % short-sighted backward pass
        ffx = X(:,t);
        timeWindow = t:min([dim.n_t-1,t]);
        dFdX = 1;
        Dy = 0;
        SiX = 0;
        for j=timeWindow
            [ffx,dF_dfx] =...
                VBA_evalFun('f',ffx,posterior.muTheta,u(:,j+1),options,dim);
            dFdX = dFdX*dF_dfx;
            [gfx,dG_dfx] =...
                VBA_evalFun('g',ffx,posterior.muPhi,u(:,j+1),options,dim);
            Gtilde = dG_dfx'*dFdX';
            Sytp1 = Qy./sigmaHat + dG_dfx'*Qi*dG_dfx./alphaHat;
            ytilde = y(:,j+1) - gfx + Gtilde*X(:,t);
            A = Gtilde'*pinv(full(Sytp1));
            SiX = SiX + A*Gtilde;
            Dy = Dy + A*(ytilde - Gtilde*m(:,t));
        end
        
        SigmaX.current{t} = VB_inv(iR{t}+SiX,indIn{t+1});
        muX(:,t) = m(:,t) + SigmaX.current{t}*Dy;
        
%         [gfx,dG_dfx] =...
%             VBA_evalFun('g',fx(:,t),posterior.muPhi,u(:,t+1),options,dim);
%         Gtilde = dG_dfx'*dF_dX{t}';
%         Sytp1 = Qy./sigmaHat + dG_dfx'*Qi*dG_dfx./alphaHat;
%         ytilde = y(:,t+1) - gfx + Gtilde*X(:,t);
%         A = Gtilde'*pinv(full(Sytp1));
%         SigmaX.current{t} = VB_inv(iR{t}+A*Gtilde,indIn{t+1});
%         muX(:,t) = m(:,t) + SigmaX.current{t}*A*(ytilde - Gtilde*m(:,t));
         

%         % prediction step: p(X_t|Y_t+1:n_t)
%         iRtilde1 = VB_inv(options.priors.SigmaX.current{t},indIn{t});
%         iRtilde2 = VB_inv(options.priors.SigmaX.current{t+1},indIn{t});
%         F = iS + (alphaHat./Te).*iQ - iRtilde2;
%         iF = VB_inv(F,indIn{t});
%         A = dF_dX{t}*(In-(alphaHat./Te).*iF);
%         iGamma = iRtilde1 - (alphaHat./Te).*(A*dF_dX{t}');
%         Gamma = VB_inv(iGamma,indIn{t});
%         muTilde(:,t) = -IN*Gamma*( ...
%             (alphaHat./Te).*A*(fx(:,t)-dF_dX{t}'*X(:,t)) ...
%             +(alphaHat./Te).*dF_dX{t}*iF*...
%             (iRtilde2*options.priors.muX(:,t+1)-iS*ohmega(:,t+1))...
%             -iRtilde1*options.priors.muX(:,t) );
%         
%         % Update step: p(X_t|Y_t:n_t)
%         iS = iGamma + C{t};
%         S = VB_inv(iS,indIn{t});
%         ohmega(:,t) = X(:,t) ...
%             + (sigmaHat./Te).*IN*S*dG_dX{t}*dy(:,t) ...
%             + S*iGamma*(muTilde(:,t)-X(:,t));
%         
%         %--- alpha-beta message (full posterior) : p(X_t|Y_1:n_t)
%         SigmaX.current{t} = VB_inv(iR{t}+iGamma-iRtilde1,indIn{t});
%         vx(:,t) = diag(SigmaX.current{t});
%         
%         muX(:,t) = m(:,t) + ...
%             IN*SigmaX.current{t}*(...
%             + iGamma*(muTilde(:,t)-m(:,t)) ...
%             - iRtilde1*(options.priors.muX(:,t)-m(:,t))...
%             );

        %-- !! only stable update: inter-time step covariance !! --%

        % Posterior inter-time covariance matrix
        A = VB_inv( (Te./alphaHat)*E -...
            (alphaHat./Te).*iQ*dF_dX{t}'*iB{t}*dF_dX{t}*iQ,indIn{t} );
        SigmaX.inter{t} = iB{t}*dF_dX{t}*iQ*A;
        
        % Entropy calculus
        jointCov = ...
            [ posterior.SigmaX.current{t+1}    posterior.SigmaX.inter{t}'
            posterior.SigmaX.inter{t}    posterior.SigmaX.current{t} ];
        indjc = [indIn{t+1}(:);indIn{t+1}(:)];
        ldj = VBA_logDet(jointCov(indjc,indjc));
        ldm = VBA_logDet(posterior.SigmaX.current{t}(indIn{t},indIn{t}));
        SX = SX + 0.5*( ldj - ldm );
        
        % update lagged covariance intermediate matrices
        iE = VB_inv(E,indIn{t+1});
        iPsi = (alphaHat./Te)*(  dF_dX{t}*iQ2*dF_dX{t}' + ...
            -(alphaHat./Te)*dF_dX{t}*iQ2*iE*iQ2*dF_dX{t}' );
        E = iPsi +(alphaHat./Te)*iQ + (sigmaHat./Te)*C{t};

        if isequal(mod(i,dim.n_t./10),0)
            set(options.display.hm(2),'string',[num2str(100*i/dim.n_t),'%']);
            drawnow
        end

    end         %--- end of backward pass ---%
end

set(options.display.hm(2),'string','OK.');
drawnow

% Gauss-Newton update step
deltaMuX = (muX - X);

% variational energy
IX = -0.5.*sigmaHat.*dy2 ...
    -0.5*alphaHat.*dx2 ...
    -0.5*sigmaHat.*Sphid2gdphi2 ...
    -0.5*alphaHat.*Sthetad2fdtheta2;
if isnan(IX) || ~isreal(IX) || div
    IX = -Inf;
end

% sufficient statistics
suffStat.SX = SX;
suffStat.Sphid2gdphi2 = Sphid2gdphi2;
suffStat.SXd2gdx2 = SXd2gdx2;
suffStat.Sthetad2fdtheta2 = Sthetad2fdtheta2;
suffStat.SXd2fdx2 = SXd2fdx2;
suffStat.SXtdfdx = SXtdfdx;
suffStat.trSx = trSx;
suffStat.gx = gx;
suffStat.vy = vy;
suffStat.dx = dx;%[muX(:,1)-fx0,muX(:,2:end) - fx];
suffStat.dy = dy;
suffStat.dx2 = dx2;
suffStat.dy2 = dy2;

% display forward and backward passes
if options.GnFigs
    try
        clf(suffStat.haf);
    catch
        suffStat.haf = figure('visible','off');
        pos = get(suffStat.haf,'position');
        set(suffStat.haf,...
            'position',pos-[pos(3)./2 1.2*pos(4) 0 0],...
            'visible','on')
    end
    h(1) = subplot(2,2,1,'parent',suffStat.haf);
    plot(h(1),mStar');
    title(h(1),'m*')
    h(2) = subplot(2,2,2,'parent',suffStat.haf);
    plot(h(2),m')
    title(h(2),'m')
    h(3) = subplot(2,2,3,'parent',suffStat.haf);
    plot(h(3),muTilde')
    title(h(3),'muTilde')
    h(4) = subplot(2,2,4,'parent',suffStat.haf);
    plot(h(4),muX')
    title(h(4),'muX')
end


