function [gx,dgdx,dgdp] = VBA_smoothNLSS(Xt,P,ut,in)
% reverses inference on (assumed AR(1)) unknown stochastic input
% function [gx,dgdx,dgdp] = VBA_smoothNLSS(Xt,P,ut,in)
% As in VBA_odeLim.m, this function evaluates the evolution/observation
% functions that are embedded in the observation function of a
% deterministic (ODE) state-space model. However, this is done here to
% enforce state noise to be AR(1), i.e. smooth.
% IN:
%   - Xt: hidden states (effectively, AR(1) state noise)
%   - P: the evolution/observation parameters
%   - ut: the input at time t
%   - in: the optional structure
% OUT:
%   - gx: the observation function evaluated at P
%   - dgdx: the derivative wrt states
%   - dgdp: the derivatives wrt the parameters
%   - d2gdxdp: the double derivatives
%------------------------------------------------------------
% Copyright (C) 2012 Jean Daunizeau / License GNU GPL v2
%------------------------------------------------------------

persistent xt dxdTheta dxdx0

% extract parameters and optional input for observation/evolution functions
if in.old.dim.n_theta > 0
    Phi = P(1:in.old.dim.n_phi);
    Theta = P(in.old.dim.n_phi+1:in.old.dim.n_phi+in.old.dim.n_theta);
else
    Phi = P(1:in.old.dim.n_phi);
    Theta = [];
end
options = in.old.options;
dim = in.old.dim;

% Check whether the system is at initial state
if isempty(xt)  % (t=0)
    dxdTheta = zeros(in.old.dim.n_theta,in.old.dim.n);
    if options.updateX0
        xt = P(in.old.dim.n_phi+in.old.dim.n_theta+1:end);
    else
        xt = in.x0;
    end
end

% apply ODE forward step:
[xt,dF_dX,dF_dP] = ...
    VBA_evalFun('f',xt,Theta,ut,options,dim);
% Add AR(1) hidden state:
xt = xt + Xt;
% apply observation mapping:
[gx,dG_dX,dG_dP] = ...
    VBA_evalFun('g',xt,Phi,ut,options,dim);


% Obtain derivatives of path wrt parameters...
dxdTheta = dF_dP + dxdTheta*dF_dX;
% ... and initial conditions
if options.updateX0
    if ~isempty(dxdx0) 
        dxdx0 = dxdx0*dF_dX;
    else % initial condition (t=0)
        dxdx0 = dF_dX;
    end
end

% Obtain derivatives of observations wrt...
dgdp = zeros(in.dim.n_phi,in.dim.p);
% ... observation parameters,
if in.old.dim.n_phi > 0
    dgdp(1:in.old.dim.n_phi,:) = dG_dP;
end
% ... evolution parameters
if in.old.dim.n_theta > 0
    dgdp(in.old.dim.n_phi+1:in.old.dim.n_phi+in.old.dim.n_theta,:) = ...
        dxdTheta*dG_dX;
end
% ... and initial conditions
if options.updateX0
    dgdp(in.old.dim.n_phi+in.old.dim.n_theta+1:end,:) = ...
        dxdx0*dG_dX;
end

% Fill in last derivatives
dgdx = dG_dX;





    
