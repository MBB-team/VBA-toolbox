function [Itheta,SigmaTheta,deltaMuTheta,suffStat] = VBA_Itheta(theta,y,posterior,suffStat,dim,u,options)
% Gauss-Newton update of the evolution parameters

if options.DisplayWin % Display progress
    set(options.display.hm(1),'string',...
        'VB Gauss-Newton on evolution parameters... ');
    set(options.display.hm(2),'string','0%');
    drawnow
end

%  Look-up which evolution parameter to update
indIn = options.params2update.theta;
indInx = options.params2update.x;

% Get precision parameters
alphaHat = posterior.a_alpha./posterior.b_alpha;

% Preallocate intermediate variables
iQx = options.priors.iQx;
Q = options.priors.SigmaTheta(indIn,indIn);
iQ = VBA_inv(Q);
muTheta0 = options.priors.muTheta;
Theta = muTheta0;
Theta(indIn) = theta;
dtheta0 = muTheta0-Theta;
dx = zeros(dim.n,dim.n_t);
div = 0;

%--- Initial condition ---%

% evaluate evolution function at current mode
[fx,dF_dX,dF_dTheta] = VBA_evalFun('f',posterior.muX0,Theta,u(:,1),options,dim,1);

% check infinite precision transition pdf
iQ2 = VBA_inv(iQx{1},indInx{1},'replace');

% posterior covariance matrix terms
d2fdx2 = dF_dTheta*iQ2*dF_dTheta';

% error terms
dx(:,1) = (posterior.muX(:,1) - fx);
dx2 = dx(:,1)'*iQ2*dx(:,1);
ddxdtheta = dF_dTheta*iQ2*dx(:,1);

%--- Loop over time series ---%
for t=1:dim.n_t-1
    
    % check infinite precision transition pdf
    iQ2 = VBA_inv(iQx{t+1},indInx{t+1},'replace');
    
    % evaluate evolution function at current mode
    [fx,dF_dX,dF_dTheta] = VBA_evalFun('f',posterior.muX(:,t),Theta,u(:,t+1),options,dim,t+1);  

    % posterior covariance matrix terms
    d2fdx2 = d2fdx2 + dF_dTheta*iQ2*dF_dTheta';

    % error terms
    dx(:,t+1) = (posterior.muX(:,t+1) - fx);
    dx2 = dx2 + dx(:,t+1)'*iQ2*dx(:,t+1);
    ddxdtheta = ddxdtheta + dF_dTheta*iQ2*dx(:,t+1);
    
    % Display progress
    if options.DisplayWin && mod(t,dim.n_t./10) < 1
        set(options.display.hm(2),'string',[num2str(floor(100*t/dim.n_t)),'%']);
        drawnow
    end
    
    % Accelerate divergent update
    if VBA_isWeird ({dx2, dF_dX, dF_dTheta})
        div = 1;
        break
    end

end

if options.DisplayWin % Display progress
    set(options.display.hm(2),'string','OK');
    drawnow
end

% posterior covariance matrix
iSigmaTheta = iQ + alphaHat.*d2fdx2(indIn,indIn);
SigmaTheta = VBA_inv(iSigmaTheta);

% mode
tmp = iQ*dtheta0(indIn) + alphaHat.*ddxdtheta(indIn);
deltaMuTheta = SigmaTheta*tmp;

% variational energy
Itheta = -0.5.*dtheta0(indIn)'*iQ*dtheta0(indIn) -0.5*alphaHat.*dx2;
if VBA_isWeird ({Itheta, SigmaTheta}) || div
    Itheta = -Inf;
end

% update sufficient statistics
suffStat.Itheta = Itheta;
suffStat.dx = dx;
suffStat.dx2 = dx2;
suffStat.dtheta = dtheta0;
suffStat.div = div;


