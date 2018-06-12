function [gx,dgdx,dgdp] = g_HRF3(Xt,P,ut,in)
% T2* contrast observation function for HRF Balloon model (log space)
% function [gx,dgdx,dgdp] = g_HRF3(Xt,P,ut,in)
% This function evaluates the hemodynamic static observation equation
% function. HEMODYNAMIC STATES ARE IN LOG SPACE.

n = size(Xt,1);

% Get parameters and states indices

try; TE = in.TE; catch, TE = 0.04; end

if isfield(in,'fullDCM') && in.fullDCM
    % estimated region-specific resting oxygen extraction fractions
    % NB: if P(in.ind1) = 0, E0 = 0.34.
%     [E0,dsdp] = sigm(P(in.ind1)-0.6633,struct('beta',1,'G0',1,'INV',0));
%     E0 = E0(:);
    E0 = 1./(1+exp(-(P(in.ind1)-0.6633)));
    dsdp = E0.*(1-E0);
    % estimated region-specific ratios of intra- to extravascular
    % components of the gradient echo signal (prior mean = 1, log-normally
    % distributed scaling factor)
    epsilon   = exp(P(in.ind2));
    % hemodynamic states indices:
    n1 = in.n1;
    n2 = in.n2;
    n3 = in.n3;
    n4 = in.n4;
    try
        n5 = in.n5;
        nreg = n./5;
    catch
        n5 = [];
        nreg = n./4;
    end
    if isfield(in,'homogeneous') && in.homogeneous
        E0 = E0.*ones(nreg,1);
        dsdp = dsdp.*ones(nreg,1);
        epsilon = epsilon.*ones(nreg,1);
        in.ind1 = [in.ind1:2:2*nreg];
        in.ind2 = [in.ind2:2:2*nreg];
    end
else
    [E0,dsdp] = VBA_sigmoid(P(1)-0.6633);
    E0 = E0(:);
    try
        epsilon = exp(P(2));
        p2 = 1;
    catch
        epsilon = 1;
        p2 = 0;
    end
    % hemodynamic states indices:
    n1 = 1;
    n2 = 2;
    n3 = 3;
    n4 = 4;
    nreg = n./4;
end

%--------------------------------------------------------------------------
% coefficients in BOLD signal model,...
V0 = 4; % resting venous volume
r0 = 25; % intravascular relaxation rate
nu0 = 40.3; % frequency offset
k1 = 4.3.*nu0.*E0.*TE;
k2 = epsilon.*r0.*E0.*TE;
k3 = 1 - epsilon;

% states ...
x3 = exp(Xt(n3,:));         % blood volume v(t)
x4 = exp(Xt(n4,:));         % dHb content q(t)

% ... and derivatives
dgdx = zeros(n,nreg);
dgdp = zeros(size(P,1),nreg);

% Evaluate observation function
x4x3 = exp(log(x4)-log(x3));
gx = V0.*(k1.*(1-x4) + k2.*(1-x4x3) + k3.*(1-x3));

% Evaluate gradient wrt states and parameters
for i=1:nreg
    if isfield(in,'fullDCM') && in.fullDCM
        if isempty(n5)
            dgdx(n1(i):n4(i),i) = [...
                0; ...
                0; ...
                -k3(i).*x3(i)+k2(i).*x4x3(i); ...
                -k1(i).*x4(i)-k2(i).*x4x3(i)].*V0;
        else
            dgdx(n1(i):n5(i),i) = [...
                0; ...
                0; ...
                -k3(i).*x3(i)+k2(i).*x4x3(i); ...
                -k1(i).*x4(i)-k2(i).*x4x3(i); 0].*V0;
        end
        dgdp(in.ind1(i),i) = V0.*4.3.*nu0.*TE.*(1-x4(i)).*dsdp(i) + ...
            V0.*epsilon(i).*r0.*TE.*(1-x4x3(i)).*dsdp(i);
        dgdp(in.ind2(i),i) = V0.*k2(i).*(1-x4x3(i)) - ...
            V0.*epsilon(i).*(1-x3(i));
    else
        dgdx(n1(i):n4(i),i) = [0; 0; ...
            -k3(i).*x3(i)+k2(i).*x4x3(i); ...
            -k1(i).*x4(i)-k2(i).*x4x3(i)].*V0;
        dgdp = V0.*4.3.*nu0.*TE.*(1-x4(i)).*dsdp(i) + ...
            V0.*epsilon(i).*r0.*TE.*(1-x4x3(i)).*dsdp(i);
        if p2
            dgdp(2) = V0.*k2(i).*(1-x4x3(i)) - ...
                V0.*epsilon(i).*(1-x3(i));
            dgdp = dgdp(:);
        end
    end
end

if isfield(in,'homogeneous') && in.homogeneous
    dgdp0 = dgdp;
    dgdp = zeros(2,nreg);
    dgdp(1,:) = sum(dgdp0(in.ind1,:),1);
    dgdp(2,:) = sum(dgdp0(in.ind2,:),1);
end

if isfield(in,'confounds') && ~isempty(in.confounds.X0)
    tsample = ut(in.confounds.indt);
    X0 = kron(eye(nreg),in.confounds.X0(tsample,:));
    gx = gx + X0*P(in.confounds.indp);
    dgdp(in.confounds.indp,:) = X0';
end


