function [IX,SigmaX,deltaMuX,suffStat] = VBA_IX_binomial(...
    X,y,posterior,suffStat,dim,u,options)
% Laplace-Kalman smoother (binomial data)

%  Look-up which hidden states to update
indIn = options.params2update.x;

% Get precision parameter
alphaHat = posterior.a_alpha./posterior.b_alpha;
iQx = options.priors.iQx;

% Preallocate intermediate variables
In = eye(dim.n);
mStar = zeros(dim.n,dim.n_t);
m = zeros(dim.n,dim.n_t);
muX = zeros(dim.n,dim.n_t);
dF_dX = cell(dim.n_t,1);
dG_dX = cell(dim.n_t,1);
iB = cell(dim.n_t,1);
iR = cell(dim.n_t,1);
iRp = cell(dim.n_t,1);
C = cell(dim.n_t,1);
SigmaX.current = cell(dim.n_t,1);
SigmaX.inter = cell(dim.n_t-1,1);
dy = zeros(dim.p,dim.n_t);
vy = zeros(dim.p,dim.n_t);
gx = zeros(dim.p,dim.n_t);
dx = zeros(dim.n,dim.n_t);
fx = zeros(dim.n,dim.n_t-1);
div = 0;

%--------------------- Forward pass -----------------------%
if options.DisplayWin
    set(options.display.hm(1),'string','VB Gauss-Newton EKF: forward pass... ');
    set(options.display.hm(2),'string','0%');
    drawnow
end

%---- Initial condition ----%
% !! Attn: MF approx. on initial condition !!

% evaluate evolution function at current mode
if ~options.priors.AR
    [fx0,dF_dX0,dF_dTheta0] = VBA_evalFun('f',posterior.muX0,posterior.muTheta,u(:,1),options,dim,1);
else
    clear VBA_smoothNLSS
    [fx0,dF_dX0,dF_dTheta0] = VBA_evalFun('f',zeros(dim.n,1),posterior.muTheta,u(:,1),options,dim,1);
end

% evaluate observation function at current mode
[gx(:,1),dG_dX{1}] = VBA_evalFun('g',X(:,1),posterior.muPhi,u(:,1),options,dim,1);

% fix numerical instabilities
gx(:,1) = checkGX_binomial(gx(:,1));

% predicted variance over binomial data
vy(:,1) = gx(:,1).*(1-gx(:,1));

% likelihood Hessian
tmp = y(:,1)./gx(:,1).^2 - (y(:,1)-1)./(1-gx(:,1)).^2;
d2gdx2 = dG_dX{1}*diag(tmp)*dG_dX{1}';

% check infinite precision transition pdf
iQ = VB_inv(iQx{1},indIn{1},'replace');
Q = VB_inv(iQ,indIn{1});
IN = diag(~~diag(iQ));

% mean-field terms
Sphid2gdphi2 = 0;
Sthetad2fdtheta2 = 0;
SXd2gdx2 = trace(dG_dX{1}*dG_dX{1}'*posterior.SigmaX.current{1});
SXd2fdx2 = trace(dF_dX0*iQ*dF_dX0'*posterior.SigmaX0);
trSx = trace(iQ*posterior.SigmaX.current{1});
SXtdfdx = 0;

% error terms
dx(:,1) = X(:,1) - fx0;
dx2 = dx(:,1)'*iQ*dx(:,1);
dy(:,1) = y(:,1) - gx(:,1);
logL = y(:,1)'*log(gx(:,1)) + (1-y(:,1))'*log(1-gx(:,1));

% covariance matrices
Rp = dF_dX0'*posterior.SigmaX0*dF_dX0 + (1./alphaHat)*Q;
C{1} = d2gdx2;
iRp{1} = VB_inv(Rp,indIn{1});
iR{1} =  iRp{1} + C{1};
R = VB_inv(iR{1},indIn{1});

% prediction : p(X_1|X_0)
mStar(:,1) = fx0;

% update : p(X_1|Y_1)
m(:,1) = mStar(:,1) + IN*R*dG_dX{1}*( dy(:,1) + dG_dX{1}'*(X(:,1)-mStar(:,1)) )./vy(:,1);

% Transient posterior
SigmaX.current{1} = R;
muX(:,1) = m(:,1);

% Entropy calculus
SX = 0.5*length(indIn{1})*log(2*pi*exp(1)) + 0.5*VBA_logDet(posterior.SigmaX.current{1},indIn{1});


%---- Sequential message-passing algorithm: forward pass (filter) ----%
for t = 2:dim.n_t
   
    % check infinite precision transition pdf
    iQ = VB_inv(iQx{t},indIn{t},'replace');
    Q =  VB_inv(iQx{t},indIn{t});
    IN = diag(~~diag(iQ));
    
    % evaluate evolution function at current mode
    [fx(:,t-1),dF_dX{t-1}] = VBA_evalFun('f',X(:,t-1),posterior.muTheta,u(:,t),options,dim,t);

    % evaluate observation function at current mode
    [gx(:,t),dG_dX{t}] = VBA_evalFun('g',X(:,t),posterior.muPhi,u(:,t),options,dim,t);
    
    % fix numerical instabilities
    gx(:,t) = checkGX_binomial(gx(:,t));
    
    % predicted variance over binomial data
    vy(:,t) = gx(:,t).*(1-gx(:,t));
    
    if isbinary(y(:,t))
        
        % likelihood Hessian
        tmp = y(:,t)./gx(:,t).^2 - (y(:,t)-1)./(1-gx(:,t)).^2;
        d2gdx2 = dG_dX{t}*diag(tmp)*dG_dX{t}';
        
        % accumulate likelihood
        logL = logL + y(:,t)'*log(gx(:,t)) + (1-y(:,t))'*log(1-gx(:,t));
        
    end
    
    % mean-field terms
    SXd2gdx2 = SXd2gdx2 + trace(dG_dX{t}*dG_dX{t}'*posterior.SigmaX.current{t});
    SXd2fdx2 = SXd2fdx2 + trace(dF_dX{t-1}*iQ*dF_dX{t-1}'*posterior.SigmaX.current{t});
    SXtdfdx = SXtdfdx - 2*trace( iQ*dF_dX{t-1}'*posterior.SigmaX.inter{t-1}' );
    trSx = trSx + trace(iQ*posterior.SigmaX.current{t});

    % error terms
    dx(:,t) = (X(:,t) - fx(:,t-1));
    dx2 = dx2 + dx(:,t)'*iQ*dx(:,t);
    if isbinary(y(:,t))
        dy(:,t) = y(:,t) - gx(:,t);
    else
        dy(:,t) = zeros(size(gx(:,t)));
    end

    % covariance matrices
    if ~options.noSXi
        B = iR{t-1} + alphaHat*dF_dX{t-1}*iQ*dF_dX{t-1}';
        iB{t-1} = VB_inv(B,indIn{t-1});
    end
    Rp = dF_dX{t-1}'*R*dF_dX{t-1} + (1./alphaHat)*Q;
    iRp{t} = VB_inv(Rp,indIn{t});
    C{t} = d2gdx2;
    iR{t} = iRp{t} + C{t};
    R = VB_inv(iR{t},indIn{t});
    
    % prediction : p(X_t|Y_1:t-1)
    mStar(:,t) = fx(:,t-1) + IN*dF_dX{t-1}'*(m(:,t-1) -X(:,t-1));     
    
    % update : p(X_t|Y_1:t)
    m(:,t) = mStar(:,t) + IN*R*dG_dX{t}*( dy(:,t) + dG_dX{t}'*(X(:,t)-mStar(:,t)) )./vy(:,1);
    
    % Transient posterior
    SigmaX.current{t} = R;
    muX(:,t) = m(:,t);


    % Display progress
    if options.DisplayWin && mod(t,dim.n_t./10) < 1
        set(options.display.hm(2),'string',[num2str(floor(100*t/dim.n_t)),'%']);
        drawnow
    end
    
    % Accelerate divergent update
    if isweird({dx2,dG_dX{t},dF_dX{t-1},SigmaX.current{t}})
        div = 1;
        break
    end
    
end         %--- end of forward pass ---%




%------------------------ backward pass ------------------------%
if options.DisplayWin  % Display progress
    set(options.display.hm(2),'string','OK');
    drawnow
    set(options.display.hm(1),'string','VB Gauss-Newton EKF: backward pass... ');
    set(options.display.hm(2),'string','0%');
    drawnow
end


%---- Sequential message-passing algorithm: backward pass (smoother) ----%
if ~div % only if forward-pass converged
    
    % Final condition beta-message (by convention)
    if ~options.noSXi
        iQ = VB_inv(iQx{dim.n_t},indIn{dim.n_t},'replace');
        E = alphaHat*iQ + sigmaHat*C{dim.n_t};
    end
    
    for i = 1:dim.n_t-1

        % time index variable change
        t = dim.n_t-i;
        
        % Posterior inter-time covariance matrix
        if ~options.noSXi
            % check infinite precision transition pdf
            iQ = VB_inv(iQx{t},indIn{t},'replace');
            iQ2 = VB_inv(iQx{t+1},indIn{t+1},'replace');
            A = VB_inv( (1./alphaHat)*E - alphaHat.*iQ*dF_dX{t}'*iB{t}*dF_dX{t}*iQ,indIn{t} );
            %         A = VB_inv( (1./alphaHat)*VB_inv(SigmaX.current{t+1},indIn{t+1})...
            %             + iQ -alphaHat.*iQ*dF_dX{t}'*iB{t}*dF_dX{t}*iQ,indIn{t} );
            %         SigmaX.inter{t} = SigmaX.inter{t}';
            %         [SigmaX] = laggedCovX(SigmaX,muX,t,...
            %             dim,iR,indIn,iRp,sigmaHat,dG_dX,dF_dX,iQy,1);
            SigmaX.inter{t} = iB{t}*dF_dX{t}*iQ*A;
            % update lagged covariance intermediate matrices
            iE = VB_inv(E,indIn{t+1});
            iPsi = alphaHat*(  dF_dX{t}*iQ2*dF_dX{t}' + -alphaHat*dF_dX{t}*iQ2*iE*iQ2*dF_dX{t}' );
            E = iPsi + alphaHat*iQ + sigmaHat*C{t};
        else
            SigmaX.inter{t} = zeros(dim.n,dim.n);
        end
        
        % Entropy calculus
        jointCov = ...
            [ posterior.SigmaX.current{t+1}  posterior.SigmaX.inter{t}
              posterior.SigmaX.inter{t}'     posterior.SigmaX.current{t} ];
        indjc = [indIn{t+1}(:);indIn{t}(:)+dim.n];
        ldj = VBA_logDet(jointCov(indjc,indjc));
        ldm = VBA_logDet(posterior.SigmaX.current{t}(indIn{t},indIn{t}));
        SX = SX + 0.5*( ldj - ldm ) + 0.5*length(indIn{t})*log(2*pi*exp(1));
        
        % display progress
        if options.DisplayWin && isequal(mod(i,dim.n_t./10),0)
            set(options.display.hm(2),'string',[num2str(100*i/dim.n_t),'%']);
            drawnow
        end
        
        % Accelerate divergent update
        if isweird(muX(:,t)) || ...
                isweird(SigmaX.current{t}) ||...
                isweird(SigmaX.inter{t})
            break
        end
        

    end         %--- end of backward pass ---%
end

if options.DisplayWin
    set(options.display.hm(2),'string','OK.');
    drawnow
end

% Gauss-Newton update step
deltaMuX = (muX - X);

% variational energy
IX = logL -0.5*alphaHat.*dx2;
if isweird({IX,SigmaX.current,SigmaX.inter}) || div
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
suffStat.dx = dx;
suffStat.dy = dy;
suffStat.dx2 = dx2;
suffStat.logL = logL;
suffStat.div = div;

% display forward and backward passes
if options.GnFigs
    try % use curretn window
        clf(suffStat.haf);
    catch % open new window
        suffStat.haf = figure('visible','off','color',[1,1,1]);
        pos = get(suffStat.haf,'position');
        set(suffStat.haf,'position',pos-[pos(3)./2 1.2*pos(4) 0 0],'visible','on')
    end
    h(1) = subplot(2,2,1,'parent',suffStat.haf);
    plot(h(1),mStar');
    title(h(1),'prediction (forward pass)')
    h(2) = subplot(2,2,2,'parent',suffStat.haf);
    plot(h(2),m')
    title(h(2),'update (forward pass)')
    h(3) = subplot(2,2,3,'parent',suffStat.haf);
    plot(h(3),dx')
    title(h(3),'state noise')
    h(4) = subplot(2,2,4,'parent',suffStat.haf);
    plot(h(4),muX')
    title(h(4),'posterior mean (backward pass)')
    axis(h,'tight')
end


