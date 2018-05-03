function [IX,SigmaX,deltaMuX,suffStat] = VBA_IX_lagged_binomial(X,y,posterior,suffStat,dim,u,options)
% lagged Laplace-EKF (Gauss-Newton update of hidden states)


%  Look-up which hidden states to update
indIn = options.params2update.x;

% Get precision parameters
alphaHat = posterior.a_alpha./posterior.b_alpha;
iQx = options.priors.iQx;

% Preallocate intermediate variables
muX = zeros(dim.n,dim.n_t);
dF_dX = cell(dim.n_t,1);
dG_dX = cell(dim.n_t,1);
dG_dPhi = cell(dim.n_t,1);
SigmaX.current = cell(dim.n_t,1);
SigmaX.inter = cell(dim.n_t-1,1);
dy = zeros(dim.p,dim.n_t);
vy = zeros(dim.p,dim.n_t);
gx = zeros(dim.p,dim.n_t);
dx = zeros(dim.n,dim.n_t);
fx = zeros(dim.n,dim.n_t-1);
div = 0;

if options.DisplayWin
    set(options.display.hm(1),'string','VB Gauss-Newton EKF: lagged forward pass... ');
    set(options.display.hm(2),'string','0%');
    drawnow
end

%---- Initial condition ----%

% form dummy posterior p(X_:0|y_:0) with augmented (lagged) states:
lag = options.backwardLag + 1;
% indx0 = dim.n*(lag-1)+1:dim.n*lag;
% m0 = zeros(dim.n*lag,1);
% m0(indx0) = posterior.muX0;
% S0 = 1e8*eye(dim.n*lag,dim.n*lag);
% S0(indx0,indx0) = posterior.SigmaX0;
m0 = repmat(posterior.muX0,lag,1);
S0 = kron(eye(lag),posterior.SigmaX0);


% evaluate evolution function at current mode
[fx0,dF_dX0] = VBA_evalFun('f',posterior.muX0,posterior.muTheta,u(:,1),options,dim,1);

% evaluate observation function at current mode
[gx(:,1),dG_dX{1},dG_dPhi{1}] = VBA_evalFun('g',X(:,1),posterior.muPhi,u(:,1),options,dim,1);

% fix numerical instabilities
gx(:,1) = VBA_finiteBinomial (gx(:,1));

% check infinite precision transition pdf
iQ = VBA_inv(iQx{1},indIn{1},'replace');

% predicted variance over binomial data
vy(:,1) = gx(:,1).*(1-gx(:,1));

% remove irregular trials
yin = find(~options.isYout(:,1));

% accumulate log-likelihood
logL = y(yin,1)'*log(gx(yin,1)) + (1-y(yin,1))'*log(1-gx(yin,1));


% error terms
dx(:,1) = X(:,1) - fx0;
dx2 = dx(:,1)'*iQ*dx(:,1);
dy(yin,1) = y(yin,1) - gx(yin,1);

% covariance matrices
GC = dG_dX{1}(:,yin)'*options.lagOp.C;
FD = dF_dX0'*options.lagOp.D;
FDC = FD - options.lagOp.C;
EuSEu = options.lagOp.Eu*S0*options.lagOp.Eu';
iEuSEu = VBA_inv(EuSEu);
EiEuSEuEu =  options.lagOp.E'*iEuSEu*options.lagOp.Eu;
EiEuSEuE = options.lagOp.E'*iEuSEu*options.lagOp.E;
xi2 = y(yin,1)./gx(yin,1).^2 - (y(yin,1)-1)./(1-gx(yin,1)).^2;
iSt = GC'*diag(xi2)*GC + alphaHat*FDC'*iQ*FDC + EiEuSEuE;
St = VBA_inv(iSt);
e1 = FD*m0 - fx0;
xi1 = dy(yin,1)./vy(yin,1);
mt = St*( GC'*(diag(xi2)*dG_dX{1}(:,yin)'*posterior.muX0-xi1) + alphaHat*FDC'*iQ*e1 + EiEuSEuEu*m0 );

%---- Sequential message-passing algorithm: lagged forward pass ----%
for t = 2:dim.n_t
    
    % check infinite precision transition pdf
    iQ = VBA_inv(iQx{t},indIn{t},'replace');
    
    % evaluate evolution function at current mode
    [fx(:,t-1),dF_dX{t-1}] = VBA_evalFun('f',X(:,t-1),posterior.muTheta,u(:,t),options,dim,t);
    
    % evaluate observation function at current mode
    [gx(:,t),dG_dX{t},dG_dPhi{t}] = VBA_evalFun('g',X(:,t),posterior.muPhi,u(:,t),options,dim,t);
    
    % fix numerical instabilities
    gx(:,t) = VBA_finiteBinomial (gx(:,t));
    
    % remove irregular trials
    yin = find(~options.isYout(:,t));
    
    % accumulate log-likelihood
    logL = logL + y(yin,t)'*log(gx(yin,t)) + (1-y(yin,t))'*log(1-gx(yin,t));
    
    % error terms
    dx(:,t) = (X(:,t) - fx(:,t-1));
    dx2 = dx2 + dx(:,t)'*iQ*dx(:,t);
    dy(yin,t) = y(yin,t) - gx(yin,t);
    
    % predicted variance over binomial data
    vy(:,t) = gx(:,t).*(1-gx(:,t));

    % covariance matrices
    GC = dG_dX{t}(:,yin)'*options.lagOp.C;
    FD = dF_dX{t-1}'*options.lagOp.D;
    FDC = FD - options.lagOp.C;
    EuSEu = options.lagOp.Eu*St*options.lagOp.Eu';
    iEuSEu = VBA_inv(EuSEu);
    EiEuSEuEu =  options.lagOp.E'*iEuSEu*options.lagOp.Eu;
    EiEuSEuE = options.lagOp.E'*iEuSEu*options.lagOp.E;
    xi2 = y(yin,t)./gx(yin,t).^2 - (y(yin,t)-1)./(1-gx(yin,t)).^2;
    iSt = GC'*diag(xi2)*GC + alphaHat*FDC'*iQ*FDC + EiEuSEuE;
    St = VBA_inv(iSt);
    e1 = dF_dX{t-1}'*X(:,t-1) - fx(:,t-1);
    xi1 = dy(yin,t)./vy(yin,t);
    mt = St*( GC'*(diag(xi2)*dG_dX{t}(:,yin)'*X(:,t)+xi1) + alphaHat*FDC'*iQ*e1 + EiEuSEuEu*mt );
    
    if t >= lag
        
        % update lagged posterior on states
        SigmaX.current{t-lag+1} = options.lagOp.M*St*options.lagOp.M';
        muX(:,t-lag+1) = options.lagOp.M*mt;
        SigmaX.inter{t-lag+1} = St(1:dim.n,dim.n+1:2*dim.n);
        
%         % Predictive density (data space)
%         V = (1./sigmaHat).*VBA_inv(iQy{t-lag+1},[]) + dG_dX{t-lag+1}'*SigmaX.current{t-lag+1}*dG_dX{t-lag+1};
%         if dim.n_phi > 0
%             V = V + dG_dPhi{t-lag+1}'*posterior.SigmaPhi*dG_dPhi{t-lag+1};
%         end
%         vy(:,t-lag+1) = diag(V);
        
    end
    
    % Display progress
    if options.DisplayWin && mod(t,dim.n_t./10) < 1
        set(options.display.hm(2),'string',[num2str(floor(100*t/dim.n_t)),'%']);
        drawnow
    end
    
    % Accelerate divergent update
    if VBA_isWeird ({dx2, dG_dX{t}, dF_dX{t-1}, SigmaX.current{t}})
        div = 1;
        break
    end
    
end         %--- end of lagged forward pass ---%



%---- End boundary of the time series ----%
for k = 2:lag
    
    % update lagged posterior on states
    ik = (k-1)*dim.n+1:k*dim.n;
    SigmaX.current{dim.n_t-(lag-k)} = St(ik,ik);
    muX(:,dim.n_t-(lag-k)) = mt(ik);
    
%     % Predictive density (data space)
%     V = (1./sigmaHat).*VBA_inv(iQy{dim.n_t-(lag-k)},[]) + dG_dX{dim.n_t-(lag-k)}'*SigmaX.current{dim.n_t-(lag-k)}*dG_dX{dim.n_t-(lag-k)};
%     if dim.n_phi > 0
%         V = V + dG_dPhi{dim.n_t-(lag-k)}'*posterior.SigmaPhi*dG_dPhi{dim.n_t-(lag-k)};
%     end
%     vy(:,dim.n_t-(lag-k)) = diag(V);
    
    if k < lag
        SigmaX.inter{dim.n_t-(lag-k)} = St(ik,ik+dim.n);
    end
    
end


if options.DisplayWin
    set(options.display.hm(2),'string','OK.');
    drawnow
end

% Gauss-Newton update step
deltaMuX = muX - X;

% variational energy
IX = logL -0.5*alphaHat.*dx2;
if VBA_isWeird ({IX, SigmaX.current, SigmaX.inter}) || div
    IX = -Inf;
end

% sufficient statistics
suffStat.IX = IX;
suffStat.gx = gx;
suffStat.vy = vy;
suffStat.dx = dx;
suffStat.dy = dy;
suffStat.logL = logL;
suffStat.dx2 = dx2;
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
    h(1) = subplot(2,1,1,'parent',suffStat.haf);
    plot(h(1),dx')
    VBA_title(h(1),'state noise')
    h(2) = subplot(2,1,2,'parent',suffStat.haf);
    plot(h(2),muX')
    VBA_title(h(2),'posterior mean')
    axis(h,'tight')
end


