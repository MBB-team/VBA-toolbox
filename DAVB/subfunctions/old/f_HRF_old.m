function [fx,dfdx,dfdp] = f_HRF_old(Xt,P,ut,in)
% Ballon (HRF) model evolution function (DCM for fMRI)
% function [fx,dfdx] = f_HRF(Xt,P,ut,in)
%
% hemodynamic response evolution function


deltat = in.deltat;
n = size(Xt,1);

% Get parameters
[E0,V0,tau0,kaf,kas,eps,a1,a2,alpha] = BOLD_parameters;
if isfield(in,'fullDCM') && in.fullDCM
    E0 = E0.*exp(P(1));
    tau0 = tau0.*exp(P(2));
    kaf = kaf.*exp(P(3));
    kas = kas.*exp(P(4));
    alpha = alpha.*exp(P(5));
    eps = 1;
else
    E0 = E0.*exp(P(1));
    tau0 = tau0.*exp(P(2));
    kaf = kaf.*exp(P(3));
    kas = kas.*exp(P(4));
    eps = eps.*exp(P(5));
    alpha = alpha.*exp(P(6));
end

% Initialize states ...
x1 = Xt(1,:);
x2 = Xt(2,:);
x3 = Xt(3,:);
x4 = Xt(4,:);
% ... and flow field, derivatives, etc...
f = zeros(n,1);
J = zeros(n,n);
dfdp = zeros(size(P,1),n);

% Evaluate flow field
f(n1) = (eps*ut - kas*x1 - kaf*(x2 - 1));
f(n2) = x1;
f(n3) = (x2 - x3.^(1./alpha))./(tau0);
f(n4) = (x2.*(1-(1-E0).^(1./x2))./E0 - x3.^(1./alpha -1).*x4)./(tau0);

% Evaluate jacobian
J(n1,:) = [ -kas , 1, 0, 0 ];
J(n2,:) = [ -kaf, 0, 1./tau0, ...
    (1 + (1-E0).^(1./x2).*((log(1-E0)./x2)-1))./(E0.*tau0)];
J(n3,:) = [ 0, 0, -(1./alpha).*x3.^(1./alpha - 1)./tau0, ...
    -(1./alpha - 1).*x4.*x3.^(1./alpha - 2)./tau0];
J(n4,:) = [ 0, 0, 0, - x3.^(1./alpha -1)./(tau0) ];

% Evaluate derivative w.r.t. parameters
dfdp(1,:) = E0.*[0,0,0,...
    ((1-E0).^(1./x2)./(1-E0)./E0-x2.*(1-(1-E0).^(1./x2))./E0.^2)./tau0];
dfdp(2,:) = [0,0,...
    -(x2-x3.^(1./alpha))./tau0,...
    -(x2.*(1-(1-E0).^(1./x2))./E0-x3.^(1./alpha-1).*x4)./tau0];
dfdp(3,:) = kaf.*[-x2+1,0,0,0];
dfdp(4,:) = kas.*[-x1,0,0,0];
if isfield(in,'fullDCM') && in.fullDCM
    dfdp(5,:) = [0,0,...
        x3.^(1./alpha)./alpha.*log(x3)./tau0,...
        x3.^(1./alpha-1)./alpha.*log(x3).*x4./tau0];
else
    dfdp(5,:) = eps.*[ut,0,0,0];
    dfdp(6,:) = [0,0,...
        x3.^(1./alpha)./alpha.*log(x3)./tau0,...
        x3.^(1./alpha-1)./alpha.*log(x3).*x4./tau0];
end

% Aply Euler discretization
fx = Xt + deltat.*f;
dfdx = eye(n) + deltat.*J;
dfdp = deltat.*dfdp;

