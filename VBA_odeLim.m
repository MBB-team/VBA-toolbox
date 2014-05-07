function [gx,out,dgdp] = VBA_odeLim(Xt,P,ut,in)
% collapses evolution and observation functions for ODE limit inversion
% function [gx,out,dgdp] = VBA_odeLim(Xt,P,ut,in)
% This function evaluates the evolution/observation functions that are
% embedded in the observation function of a deterministic (ODE) state-space
% model. It also outputs the derivatives with respect to the parameters.
% This is an overloaded function, since it is used during the
% initialization of the full stochastic state-space model inversion.
% IN:
%   - Xt: [irrelevant]
%   - P: the evolution/observation parameters
%   - ut: the input at time t
%   - in: the optional structure
% OUT:
%   - gx: the observation function evaluated at P
%   - out: [used to pass gradients of hidden states dynamics wrt evolution parameter]
%   - dgdp: the derivatives wrt the parameters

persistent t xt dxdTheta dxdx0

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
if isempty(t)
    t = 1;
    dxdTheta = zeros(dim.n_theta,dim.n);
    if options.updateX0
        xt = P(dim.n_phi+dim.n_theta+1:end);
    else
        xt = options.priors.muX0;
    end
else
    t = t+1;
end

% apply ODE forward step:
[xt,dF_dX,dF_dP] = VBA_evalFun('f',xt,Theta,ut,options,dim,t);
[gx,dG_dX,dG_dP] = VBA_evalFun('g',xt,Phi,ut,options,dim,t);

% Obtain derivatives of path wrt parameters...
if dim.n_theta > 0
    dxdTheta = dF_dP + dxdTheta*dF_dX;
end
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
if dim.n_phi > 0
    dgdp(1:dim.n_phi,:) = dG_dP;
end
% ... evolution parameters
if dim.n_theta > 0
    dgdp(dim.n_phi+1:dim.n_phi+dim.n_theta,:) = dxdTheta*dG_dX;
end
% ... and initial conditions
if options.updateX0
    dgdp(dim.n_phi+dim.n_theta+1:end,:) = dxdx0*dG_dX;
end


% for hidden states book keeping
out.xt = xt;
out.dx = [zeros(dim.n_phi,dim.n);dxdTheta;dxdx0];
out.t = t;


