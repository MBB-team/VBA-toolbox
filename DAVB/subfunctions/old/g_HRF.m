function [gx,dgdx,dgdp] = g_HRF(Xt,P,ut,in)
% T2* contrast observation function for HRF Balloon model (DCM for fMRI)
% function [gx,dgdx,dgdp] = g_HRF(Xt,P,ut,in)
% This function evaluates the hemodynamic static observation equation
% function.

n = size(Xt,1);


% Get parameters
[E0,V0,tau0,kaf,kas,eps,alpha,k1,k2,k3] = BOLD_parameters;
if isfield(in,'fullDCM') && in.fullDCM
    V0 = V0.*exp(P(in.ind2));
    k1 = 7*E0.*exp(P(in.ind1));
    k3 = 2*E0.*exp(P(in.ind1)) - 0.2;
    n1 = in.n1;
    n2 = in.n2;
    n3 = in.n3;
    n4 = in.n4;
    n5 = in.n5;
    nreg = n./5;
else
    V0 = V0.*exp(P(1));
    n1 = 1;
    n2 = 2;
    n3 = 3;
    n4 = 4;
    nreg = n./4;
end
k2 = 2;

% Initialize states ...
x3 = Xt(n3,:);         % blood volume v(t)
x4 = Xt(n4,:);         % dHb content q(t)
% ... and derivatives
dgdx = sparse(n,nreg);
dgdp = sparse(size(P,1),nreg);

% Evaluate observation function
x4x3 = exp(log(x4)-log(x3));
gx = V0.*(k1.*(1-x4) + k2.*(1-x4x3) + k3.*(1-x3));


% Evaluate gradient wrt states and parameters
for i=1:nreg
    if isfield(in,'fullDCM') && in.fullDCM
        dgdx(n1(i):n5(i),i) = ...
            [0; 0; -k3(i)+k2.*x4(i)./x3(i).^2; -k1(i)-k2./x3(i); 0].*V0(i);
        dgdp(in.ind1(i),i) = V0(i).*(7*E0.*exp(P(in.ind1(i))).*(1-x4(i))...
            +2*E0.*exp(P(in.ind1(i))).*(1-x3(i)));
        dgdp(in.ind2(i),i) = gx(i);
    else
        dgdx(n1(i):n4(i),i) = ...
            [0; 0; -k3(i)+k2.*x4(i)./x3(i).^2; -k1(i)-k2./x3(i)].*V0(i);
        dgdp = gx;
    end
end



