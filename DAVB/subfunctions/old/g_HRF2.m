function [gx,dgdx,dgdp] = g_HRF2(Xt,P,ut,in)
% T2* contrast observation function for HRF Balloon model (DCM for fMRI)
% function [gx,dgdx,dgdp] = g_HRF2(Xt,P,ut,in)
% This function evaluates the hemodynamic static observation equation
% function.

n = size(Xt,1);


try % EPI echo time
    TE = M.TE;
catch
    TE = 0.04;
end
% resting venous volume
V0        = 100*0.04;
% slope r0 of intravascular relaxation rate R_iv as a function of oxygen
% saturation Y:  R_iv = r0*[(1-Y)-(1-Y0)]
r0        = 25; % [Hz]
% frequency offset at the outer surface of magnetized vessels
nu0       = 40.3; % [Hz]


if isfield(in,'fullDCM') && in.fullDCM
    % estimated region-specific resting oxygen extraction fractions
    % NB: if P(in.ind1) = 0, E0 = 0.34.
    [E0,dsdp] = sigm(P(in.ind1)-0.6633,struct('beta',2));
%     E0        = 0.34 + P(in.ind1);
    
    % estimated region-specific ratios of intra- to extravascular components of
    % the gradient echo signal (prior mean = 1, log-normally distributed
    % scaling factor)
    epsilon   = exp(P(in.ind2));
else
    [E0,dsdp] = sigm(P(1)-0.6633,struct('beta',2));
%     E0 = 0.34 + P(1);
    try
        epsilon = exp(P(2));
        p2 = 1;
    catch
        epsilon = 1;
        p2 = 0;
    end
end

%--------------------------------------------------------------------------

% coefficients in BOLD signal model
E0 = E0(:);
k1       = 4.3.*nu0.*E0.*TE;
k2       = epsilon.*r0.*E0.*TE;
k3       = 1 - epsilon;


% Get parameters
if isfield(in,'fullDCM') && in.fullDCM
    n1 = in.n1;
    n2 = in.n2;
    n3 = in.n3;
    n4 = in.n4;
    n5 = in.n5;
    nreg = n./5;
else
    n1 = 1;
    n2 = 2;
    n3 = 3;
    n4 = 4;
    nreg = n./4;
end

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
x4x32 = exp(log(x4)-2*log(x3));
for i=1:nreg
    if isfield(in,'fullDCM') && in.fullDCM
        dgdx(n1(i):n5(i),i) = ...
            [0; 0; -k3(i)+k2(i).*x4x32(i); -k1(i)-k2(i)./x3(i); 0].*V0;
        dgdp(in.ind1(i),i) = V0.*4.3.*nu0.*TE.*(1-x4(i)).*dsdp(i) + ...
            V0.*epsilon(i).*r0.*TE.*(1-x4x3(i)).*dsdp(i);
        dgdp(in.ind2(i),i) = V0.*k2(i).*(1-x4x3(i)) - ...
            V0.*epsilon(i).*(1-x3(i));
    else
        dgdx(n1(i):n4(i),i) = ...
            [0; 0; -k3(i)+k2(i).*x4x32(i); -k1(i)-k2(i)./x3(i)].*V0;
        dgdp = V0.*4.3.*nu0.*TE.*(1-x4(i)).*dsdp(i) + ...
            V0.*epsilon(i).*r0.*TE.*(1-x4x3(i)).*dsdp(i);
        if p2
            dgdp(2) = V0.*k2(i).*(1-x4x3(i)) - ...
                V0.*epsilon(i).*(1-x3(i));
            dgdp = dgdp(:);
        end
    end
end



