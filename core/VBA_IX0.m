
function [IX0,SigmaX0,deltaMuX0,suffStat] = VBA_IX0(X0,y,posterior,suffStat,dim,u,options)
% Gauss-Newton update of initial conditions



if options.DisplayWin % Display progress
    set(options.display.hm(1),'string','VB Gauss-Newton on initial conditions... ');
    set(options.display.hm(2),'string',' ');
    drawnow
end

% check infinite precision transition pdf
iQx0 = VBA_inv(options.priors.iQx{1},options.params2update.x{1},'replace');
IN = diag(~~diag(iQx0));

% Get precision parameters
alphaHat = posterior.a_alpha./posterior.b_alpha;

% Preallocate intermediate variables
indIn = options.params2update.x0;
muX0 = options.priors.muX0;
x0 = muX0;
x0(indIn) = X0;

% Evaluate evolution function at current mode
[fx0,dF_dX0] = VBA_evalFun('f',x0,posterior.muTheta,u(:,1),options,dim,1);

% error terms
dx = IN*(posterior.muX(:,1)- fx0);
dx2 = dx'*iQx0*dx;
dx0 = muX0-x0;

% posterior covariance matrix terms
Q = options.priors.SigmaX0(indIn,indIn);
iQ = VBA_inv(Q);
iSigmaX0 = iQ + alphaHat.*dF_dX0(indIn,:)*iQx0(indIn,indIn)*dF_dX0(indIn,:)';

% posterior covariance matrix
SigmaX0 = VBA_inv(iSigmaX0);

% mode
tmp = iQ*dx0(indIn) + alphaHat.*dF_dX0(indIn,:)*iQx0(indIn,indIn)*dx;
deltaMuX0 = SigmaX0*tmp;

% variational energy
IX0 = -0.5.*dx0'*iQ*dx0 - 0.5*alphaHat.*dx2;
if VBA_isWeird (IX0)
    div = 1;
    IX0 = -Inf;
else
    div = 0;
end

% update sufficient statistics
dx20 = suffStat.dx(:,1)'*iQx0*suffStat.dx(:,1);
suffStat.dx2 = suffStat.dx2 - dx20 + dx2; % correct states squared error
suffStat.dx(:,1) = dx;
suffStat.dx0 = dx0;
suffStat.div = div;
suffStat.IX0 = IX0;


