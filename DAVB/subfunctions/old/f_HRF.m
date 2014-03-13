function [fx,dfdx,dfdp] = f_HRF(Xt,P,ut,in)
% Balloon (HRF) model evolution function
% function [fx,dfdx] = f_HRF(Xt,P,ut,in)
% This function evaluates the evolution function derived from the balloon
% model for the hemodynamic response function. It can be called in two
% ways: (i) as a "stand-alone" evolution function, whereby the system's
% states are hemodynamic states of the balloon model, or (ii) as a
% generalized observation function, where the real system's states are the
% neuronal states of a DCM for fMRI model.


deltat = in.deltat;
n = size(Xt,1);

% Get parameters
[E0,V0,tau0,kaf,kas,eps,alpha] = BOLD_parameters;
if isfield(in,'fullDCM') && in.fullDCM
    nreg = n./5;
    ind1 = in.ind1;
    ind2 = in.ind2;
    ind3 = in.ind3;
    ind4 = in.ind4;
    n1 = in.n1;
    n2 = in.n2;
    n3 = in.n3;
    n4 = in.n4;
    n5 = in.n5;
    eps = 1;
    alpha = alpha.*exp(P(in.ind5));
else
    nreg = n./4;
    ind1 = 1;
    ind2 = 2;
    ind3 = 3;
    ind4 = 4;
    n1 = 1;
    n2 = 2;
    n3 = 3;
    n4 = 4;
    eps = eps.*exp(P(5));
    alpha = alpha.*exp(P(6));
end
[E0,dsdp] = sigm(P(ind1)-0.6633,struct('beta',2));
E0 = E0(:);
% E0 = E0.*exp(P(ind1));
tau0 = tau0.*exp(P(ind2));
kaf = kaf.*exp(P(ind3));
kas = kas.*exp(P(ind4));

% Initialize states ...
x1 = Xt(n1,:);      % vasodilatory signal s(t)
x2 = Xt(n2,:);      % blood inflow f(t)
x3 = Xt(n3,:);      % blood volume v(t)
x4 = Xt(n4,:);      % dHb content q(t)
% ... and flow field, derivatives, etc...
f = zeros(n,1);
J = sparse(n,n);
dfdp = sparse(size(P,1),n);

% Evaluate flow field
f(n1) = (eps.*ut - kas.*x1 - kaf.*(x2 - 1));
f(n2) = x1;
f(n3) = (x2 - x3.^(1./alpha))./(tau0);
f(n4) = (x2.*(1-(1-E0).^(1./x2))./E0 - x3.^(1./alpha -1).*x4)./(tau0);

% Evaluate jacobian and gradients wrt parameters
for i=1:nreg
    J(n1(i),n1(i):n1(i)+3) = [ -kas(i) , 1, 0, 0 ];
    J(n2(i),n1(i):n1(i)+3) = [ -kaf(i), 0, 1./tau0(i), ...
        (1 + (1-E0(i)).^(1./x2(i)).*((log(1-E0(i))./x2(i))-1))./(E0(i).*tau0(i))];
    J(n3(i),n1(i):n1(i)+3) = [ 0, 0, -(1./alpha(i)).*x3(i).^(1./alpha(i) - 1)./tau0(i), ...
        -(1./alpha(i) - 1).*x4(i).*x3(i).^(1./alpha(i) - 2)./tau0(i)];
    J(n4(i),n1(i):n1(i)+3) = [ 0, 0, 0, - x3(i).^(1./alpha(i) -1)./(tau0(i)) ];
    
    dfdp(ind1(i),n1(i):n1(i)+3) = dsdp(i).*[0,0,0,...
        ((1-E0(i)).^(1./x2(i))./(1-E0(i))./E0(i)-x2(i).*(1-(1-E0(i)).^(1./x2(i)))./E0(i).^2)./tau0(i)];
    dfdp(ind2(i),n1(i):n1(i)+3) = [0,0,...
        -(x2(i)-x3(i).^(1./alpha(i)))./tau0(i),...
        -(x2(i).*(1-(1-E0(i)).^(1./x2(i)))./E0(i)-x3(i).^(1./alpha(i)-1).*x4(i))./tau0(i)];
    dfdp(ind3(i),n1(i):n1(i)+3) = kaf(i).*[-x2(i)+1,0,0,0];
    dfdp(ind4(i),n1(i):n1(i)+3) = kas(i).*[-x1(i),0,0,0];
    
    if isfield(in,'fullDCM') && in.fullDCM
        J(n5(i),n1(i)) = eps;
        dfdp(in.ind5(i),n1(i):n1(i)+3) = [0,0,...
            x3(i).^(1./alpha(i))./alpha(i).*log(x3(i))./tau0(i),...
            x3(i).^(1./alpha(i)-1)./alpha(i).*log(x3(i)).*x4(i)./tau0(i)];
    else
        dfdp(5,:) = eps.*[ut,0,0,0];
        dfdp(6,:) = [0,0,...
            x3.^(1./alpha)./alpha.*log(x3)./tau0,...
            x3.^(1./alpha-1)./alpha.*log(x3).*x4./tau0];
    end
    
    
end


% Apply Euler discretization
fx = Xt + deltat.*f;
dfdx = eye(n) + deltat.*J;
dfdp = deltat.*dfdp;

