function [fx,dfdx,dfdp] = f_HRF2(Xt,P,ut,in)
% Balloon (HRF) model evolution function in log-space
% function [fx,dfdx,dfdp] = f_HRF2(Xt,P,ut,in)
% This function evaluates the evolution function derived from the balloon
% model for the hemodynamic response function. It can be called in two
% ways: (i) as a "stand-alone" evolution function, whereby the system's
% states are hemodynamic states of the balloon model, or (ii) as a
% generalized observation function, where the real system's states are the
% neuronal states of a DCM for fMRI model. Note that the hemodynamic states
% are in log-space, for positivity constrints.


deltat = in.deltat;
n = size(Xt,1);

% Get parameters
[E0,V0,tau0,kaf,kas,epsilon,alpha] = BOLD_parameters;
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
    try, n5=in.n5; catch, n5=[];end
    epsilon = 1;
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
    epsilon = epsilon.*exp(P(5));
    alpha = alpha.*exp(P(6));
end
% [E0,dsdp] = sigm(P(ind1)-0.6633,struct('beta',1,'G0',1,'INV',0));
% E0 = E0(:);
E0 = 1./(1+exp(-(P(ind1)-0.6633)));
dsdp = E0.*(1-E0);
tau0 = tau0.*exp(P(ind2));
kaf = kaf.*exp(P(ind3));
kas = kas.*exp(P(ind4));

if isfield(in,'xshift') % for numerical stability
    xshift = in.xshift;
else
    xshift = 0;
end

% hemodynamic states ...
if isfield(in,'linearized') && in.linearized
    x1 = zeros(nreg,1);
    x2 = ones(nreg,1);
    x3 = ones(nreg,1);
    x4 = ones(nreg,1);
else
    % vasodilatory signal s(t)
    x1 = Xt(n1,:);
    % blood inflow f(t)
    if isfield(in,'logx2') && ~in.logx2
        x2 = Xt(n2,:) + 1; % deviation to steady-state!
    else
        x2 = exp(Xt(n2,:)) + xshift;
    end
    % blood volume v(t)
    x3 = exp(Xt(n3,:)) + 0.*xshift;
    % dHb content q(t)
    x4 = exp(Xt(n4,:)) + 0.*xshift;
end
% blood outflow
fv = x3.^(1./alpha);
% d[blood flow]/dXt(3)
dfvdx = (1./alpha).*fv;
% oxygen extraction
ff = (1-(1-E0).^(1./x2))./E0;
% d[O2 extraction]/dXt(2)
if isfield(in,'logx2') && ~in.logx2
    dffdx = log(1-E0).*(1-E0).^(1./x2)./(E0.*x2.^2);
else
    dffdx = log(1-E0).*(1-E0).^(1./x2)./(E0.*x2);
end
% ... and flow field, derivatives, etc...
f = zeros(n,1);
J = zeros(n,n);
dfdp = zeros(size(P,1),n);

% Evaluate flow field
f(n1) = (epsilon.*ut - kas.*x1 - kaf.*(x2 - 1));
if isfield(in,'logx2') && ~in.logx2
    f(n2) = x1;
else
    f(n2) = x1./x2;
end
f(n3) = (x2 - fv)./(tau0.*x3);
f(n4) = (x2.*ff./x4 - fv./x3)./tau0;


% Evaluate jacobian and gradients wrt parameters
for i=1:nreg
    
    if isfield(in,'logx2') && ~in.logx2
        J(n1(i),n1(i):n1(i)+3) = [ -kas(i) , 1, 0, 0 ];
        J(n2(i),n1(i):n1(i)+3) = [ -kaf(i), 0, ...
            1./(tau0(i).*x3(i)), ...
            (ff(i)+x2(i).*dffdx(i))./(x4(i).*tau0(i))];
    else
        J(n1(i),n1(i):n1(i)+3) = [ -kas(i) , 1./x2(i), 0, 0 ];
        J(n2(i),n1(i):n1(i)+3) = [ -kaf(i).*x2(i), -x1(i)./x2(i), ...
            x2(i)./(tau0(i).*x3(i)), ...
            x2(i).*(ff(i)+dffdx(i))./(tau0(i).*x4(i))];
    end
    J(n3(i),n1(i):n1(i)+3) = [ 0, 0,...
        -f(n3(i)) - dfvdx(i)./(tau0(i).*x3(i)) , ...
        (fv(i)-dfvdx(i))./(tau0(i).*x3(i))];
    J(n4(i),n1(i):n1(i)+3) = [ 0, 0, 0, ...
        -(x2(i).*ff(i))./(tau0(i).*x4(i)) ];
    
    
    if isfield(in,'linearized') && in.linearized
        
        tmp = log(1-E0(i))./E0(i);
        
        dfdp(ind1(i),n1(i):n1(i)+3) = [0,0,0,...
            -dsdp(i).*(1+tmp).*Xt(n2(i))./(tau0(i).*E0(i))];
        dfdp(ind2(i),n1(i):n1(i)+3) = [0,0,...
            (-Xt(n2(i))+Xt(n3(i))./alpha(i))./tau0(i),...
            -((ff(i)+dffdx(i)).*Xt(n2(i)) + (fv(i)-dfvdx(i)).*Xt(n3(i)) ...
            - ff(i).*Xt(n4(i)))./tau0(i)];
        dfdp(ind3(i),n1(i):n1(i)+3) = kaf(i).*[-Xt(n2(i)),0,0,0];
        dfdp(ind4(i),n1(i):n1(i)+3) = kas(i).*[-Xt(n1(i)),0,0,0];
        
        if isfield(in,'fullDCM') && in.fullDCM
            if  ~isempty(n5)
                J(n5(i),n1(i)) = epsilon;
            end
            dfdp(in.ind5(i),n1(i):n1(i)+3) = [0,0,...
                Xt(n3(i))./(tau0(i).*alpha(i)),...
                Xt(n3(i))./(tau0(i).*alpha(i))];
        else
            dfdp(5,:) = epsilon.*[ut,0,0,0];
            dfdp(6,:) = [0,0,...
                Xt(n3(i))./(tau0(i).*alpha(i)),...
                Xt(n3(i))./(tau0(i).*alpha(i))];
        end
        
    else
        
        % gradient wrt parameters
        dfdp(ind1(i),n1(i):n1(i)+3) = [0,0,0,...
            (((1-E0(i)).^(-1+1./x2(i)))-x2(i).*ff(i))...
            .*dsdp(i)./(tau0(i).*x4(i).*E0(i))];
        dfdp(ind2(i),n1(i):n1(i)+3) = [0,0,...
            -(x2(i) - fv(i))./(tau0(i).*x3(i)),...
            -(x2(i).*ff(i)./x4(i) - fv(i)./x3(i))./tau0(i)];
        dfdp(ind3(i),n1(i):n1(i)+3) = kaf(i).*[-x2(i)+1,0,0,0];
        dfdp(ind4(i),n1(i):n1(i)+3) = kas(i).*[-x1(i),0,0,0];
        
        % complement gradients if used for full DCM model
        if isfield(in,'fullDCM') && in.fullDCM
            if  ~isempty(n5)
                J(n5(i),n1(i)) = epsilon;
            end
            dfdp(in.ind5(i),n1(i):n1(i)+3) = [0,0,...
                log(x3(i)).*fv(i)./(tau0(i).*x3(i).*alpha(i)),...
                log(x3(i)).*fv(i)./(tau0(i).*x3(i).*alpha(i))];
        else
            dfdp(5,:) = epsilon.*[ut,0,0,0];
            dfdp(6,:) = [0,0,...
                log(x3).*fv./(tau0.*x3.*alpha),...
                log(x3).*fv./(tau0.*x3.*alpha)];
        end
        
    end
    
end

% Apply Euler discretization
if isfield(in,'linearized') && in.linearized
     if isfield(in,'fullDCM') && in.fullDCM && ~isempty(n5)
         Cu = 0;
     else         
         Cu = zeros(n,1);
         Cu(n1) = epsilon.*ut;
     end
    fx = Xt + deltat.*(J'*Xt + Cu);
else
    fx = Xt + deltat.*f;
end
dfdp = deltat.*dfdp;
dfdx = eye(n) + deltat.*J;


