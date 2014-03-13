function [gx,dgdx,dgdp] = g_HRF_old(Xt,P,ut,in)
% T2* contrast observation function for HRF Balloon model (DCM for fMRI)
% function [gx,dgdx,dgdp] = g_HRF_old(Xt,P,ut,in)
% This function evaluates the hemodynamic static observation equation
% function.

n = size(Xt,1);
nreg = n./4;

% Get parameters
[E0,V0,tau0,kaf,kas,eps,a1,a2,alpha] = BOLD_parameters;
if isfield(in,'fullDCM') && in.fullDCM
    V0 = V0.*exp(P(in.ind6));
    a1 = 7*E0.*exp(P(in.ind1)) +2;
    a2 = 2.2 -2*E0.*exp(P(in.ind1));
    n1 = in.n1;
    n2 = in.n2;
    n3 = in.n3;
    n4 = in.n4;
else
    V0 = V0.*exp(P(1));
    n1 = 1;
    n2 = 2;
    n3 = 3;
    n4 = 4;
end

% Initialize states ...
x3 = Xt(n3,:);         % q(t)
x4 = Xt(n4,:);         % v(t)
% ... and derivatives
dgdx = sparse(n,nreg);
dgdp = sparse(size(P,1),nreg);

% Evaluate obseravtion function
gx = V0.*(a1.*(1-x4) - a2.*(1-x3));

% dgdx = [];
% dgdp = [];

% Evaluate gradient wrt states and parameters
for i=1:nreg
    dgdx(n1(i):n4(i),i) = [0; 0; a2(i); -a1(i)].*V0(i);
    if isfield(in,'fullDCM') && in.fullDCM
        dgdp(in.ind1(i),i) = V0(i).*(7*E0.*exp(P(in.ind1(i))).*(1-x4(i)) ...
            +2*E0.*exp(P(in.ind1(i))).*(1-x3(i)));
        dgdp(in.ind6(i),i) = gx(i);
    else
        dgdp = gx;
    end
end



