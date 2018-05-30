function [muX,SigmaX,suffStat] = VBA_EKF(y,u,posterior,dim,options,flag)
% standard EKF & computation of the predictive density
% function [muX,SigmaX,suffStat] = VBA_EKF(y,u,posterior,dim,options,flag)
%
% This function inverts any nonlinear state-space model of the form:
%   y_t = g( x_t,u_t,phi ) + e_t
%   x_t+1 = f( x_t,u_t,theta ) + f_t
% using a standard extended Kalman filter (EKF).
%
% IN :          [ see VBA_NLStateSpaceModel.m ]
%   - y: pxn_t mesurements matrix
%   - u: mxn_t known input matrix (which is required as an argument in
%   the obvservation/evolution functions)
%        
%   - posterior: posterior pdf structure (see VBA_NLStateSpaceModel.m).
%   This is used to extract first and second order parameters of the
%   generative model (theta, phi, alpha, sigma), as well as the initial
%   conditions of the hidden states (X0).
%   - dim: a structure variable containing the dimensions of the 3 sets of
%   model unknown variables (see VBA_NLStateSpaceModel.m)
%   - options: user-defined structure containing specific informations
%   regarding the model, ie (see VBA_check.m):
%       .f_fname (resp. g_fname): name/handle of the function that outputs
%       the evolution (resp. observation) of the hidden states.
%       .u0: the mx1 initial value of the input {0}
%       .inF: a (possibly structure) variable containing the additional
%       (internal) fixed parameters which may have to be sent to the
%       evolution function (eg pointing to different variants) {[]}
%       .inG: idem for the observation function {[]}
%   -  flag: a switch for computing just the mode of the deterministic
%   predictive density (flag=0), the 1st and 2d order statistics of the EKF
%   ({flag=1}), or the 1st and 2d order statistics of the predictive
%   density (flag = 2).
%
% OUT:
%   - muX: posterior mean of the hidden states X (nxn_t matrix)
%   - SigmaX: covariance matrices of the variational posterior pdf of
%       the dynamic hidden-states.

if numel(options.sources) > 1
    error('*** EKF is not yet compatible with multisource observations');
end

if options.sources.type > 1
    error('*** EKF is not yet compatible with multinomial observations');
end

% By default, this function implements an EKF:
if ~exist('flag','var') || isempty(flag)
    flag = 1;
end

% This checks and fills in required dummy variables
if isempty(u)
    u = zeros(1,dim.n_t);
end
[dim.p,dim.n_t] = size(y);
try
    dim.u;
catch
    dim.u = size(u,1);
end
if isfield(options,'microU') && options.microU
    u = VBA_getU(u,options,dim,'2macro');
end
if ~isfield(options,'nout_f')
    options.nout_f = nargout(options.f_fname);
end
if ~isfield(options,'nout_g')
    options.nout_g = nargout(options.g_fname);
end
if ~isfield(options,'OnLine')
    options.OnLine = 0;
end
try
    X0 = posterior.muX0;
    SigmaX0 = posterior.SigmaX0;
catch
    X0 = zeros(dim.n,1);
    SigmaX0 = zeros(dim.n,dim.n);
end
try
    theta = posterior.muTheta;
catch
    theta = [];
end
try
    phi = posterior.muPhi;
catch
    phi = [];
end
try
    iQx = options.priors.iQx;
    iQy = options.priors.iQy;
catch
    iQx = cell(dim.n_t,1);
    iQy = cell(dim.n_t,1);
    for t= 1:dim.n_t
        iQx{t} = eye(dim.n);
        iQy{t} = eye(dim.p);
    end
end


switch flag
    case 0
        suffStat = [];
        str = 'deterministic time series';
    case {1,2}
        try
            alpha = posterior.a_alpha(end)./posterior.b_alpha(end);
            if options.sources.type == 0
                sigma = posterior.a_sigma(end)./posterior.b_sigma(end);
            end
        catch
            error('Not enough info in posterior structure!')
        end
        mStar = zeros(dim.n,dim.n_t);
        SigmaX = cell(dim.n_t,1);
        %--- Initialize sufficient statistics time-series ---%
        suffStat = VBA_getSuffStat(options);
        if isequal(flag,1)
            str = 'standard EKF';
        else
            str = 'predictive density';
        end
end
muX = zeros(dim.n,dim.n_t);
gx = zeros(dim.p,dim.n_t);


% First time iteration (from initial conditions)
if ~options.OnLine && options.verbose
    fprintf(1,['Deriving ',str,' ...'])
end
if flag>=1
    %--- Prediction
    [fx0,dF_dX0] = VBA_evalFun('f',X0,theta,u(:,1),options,dim,1);
    mStar(:,1) = fx0;
    Rp = dF_dX0'*SigmaX0*dF_dX0 + 1./alpha.*VBA_inv(iQx{1});
    if flag == 1 % EKF update
        [gx(:,1),dG_dX] = VBA_evalFun('g',mStar(:,1),phi,u(:,1),options,dim,1);
        iRp = pinv(Rp);
        C =  dG_dX*iQy{1}*dG_dX';
        iSX = iRp + sigma*C;
        SigmaX{1} = pinv( iSX );
        muX(:,1) = mStar(:,1) + sigma.*SigmaX{1}*dG_dX*iQy{1}* (y(:,1)-gx(:,1));
    else % Predictive density
        muX(:,1) = mStar(:,1);
        SigmaX{1} = Rp;
    end
    % get predicted observation at the mode
    [gx(:,1),dG_dX] = VBA_evalFun('g',muX(:,1),phi,u(:,1),options,dim,1);
    suffStat.dy(:,1) = y(:,1) - gx(:,1);
    if options.sources.type == 0
        suffStat.vy(:,1) = diag( sigma.^-1.*pinv(iQy{1}) + dG_dX'*SigmaX{1}*dG_dX );
        suffStat.dy2 = suffStat.dy2 + suffStat.dy(:,1)'*iQy{1}*suffStat.dy(:,1);
    else
        suffStat.vy(:,1) = gx(:,1).*(1-gx(:,1));
        suffStat.logL = y(:,1)'*log(gx(:,1)) + (1-y(:,1))'*log(1-gx(:,1));
    end
    suffStat.dx(:,1) = muX(:,1) - fx0;
    suffStat.dx2 = suffStat.dx2 + suffStat.dx(:,1)'*iQx{1}*suffStat.dx(:,1);
else
    muX(:,1) = VBA_evalFun('f',X0,theta,u(:,1),options,dim);
end


% Loop over time samples
if ~options.OnLine && options.verbose
    fprintf(1,'%6.2f %%',0)
end
for t = 1:dim.n_t-1
    if flag >= 1
        %-- Prediction
        [fx,dF_dX] = VBA_evalFun('f',muX(:,t),theta,u(:,t+1),options,dim,t+1);
        mStar(:,t+1) = fx;
        Rp = dF_dX'*SigmaX{t}*dF_dX + 1./alpha.*VBA_inv(iQx{t+1});
        if flag == 1    % EKF update
            [gx(:,t+1),dG_dX] = VBA_evalFun('g',mStar(:,t+1),phi,u(:,t+1),options,dim,t+1);
            C =  dG_dX*iQy{t+1}*dG_dX';
            iRp = pinv(Rp);
            iSX = iRp + sigma*C;
            SigmaX{t+1} = pinv( iSX );
            muX(:,t+1) = mStar(:,t+1) + sigma.*SigmaX{t+1}*dG_dX*iQy{t+1}* (y(:,t+1)-gx(:,t+1));
        else
            muX(:,t+1) = mStar(:,t+1);
            SigmaX{t+1} = Rp;
        end
        % get predicted observation at the mode
        [gx(:,t+1),dG_dX] = VBA_evalFun('g',muX(:,t+1),phi,u(:,t+1),options,dim,t+1);
        suffStat.dy(:,t+1) = y(:,t+1) - gx(:,t+1);
        if options.sources.type == 0
            suffStat.vy(:,t+1) = diag( sigma.^-1.*pinv(iQy{t+1}) + dG_dX'*SigmaX{t+1}*dG_dX );
            suffStat.dy2 = suffStat.dy2 + suffStat.dy(:,t+1)'*iQy{t+1}*suffStat.dy(:,t+1);
        else
            suffStat.vy(:,t+1) = gx(:,t+1).*(1-gx(:,t+1));
            suffStat.logL = suffStat.logL + y(:,t+1)'*log(gx(:,t+1)) + (1-y(:,t+1))'*log(1-gx(:,t+1));
        end
        suffStat.dx(:,t+1) = muX(:,t+1) - fx;
        suffStat.dx2 = suffStat.dx2 + suffStat.dx(:,t+1)'*iQx{t+1}*suffStat.dx(:,t+1);
    else
        muX(:,t+1) = VBA_evalFun('f',muX(:,t),theta,u(:,t+1),options,dim);
    end
    if ~options.OnLine && isequal(mod(t,32),0) && options.verbose
        fprintf(1,repmat('\b',1,8))
        fprintf(1,'%6.2f %%',100*t/dim.n_t)
    end
end
if ~options.OnLine  && options.verbose
    fprintf(1,repmat('\b',1,8))
    fprintf(' OK.')
    fprintf('\n')
end
if flag >= 1
    suffStat.gx = gx;
end

