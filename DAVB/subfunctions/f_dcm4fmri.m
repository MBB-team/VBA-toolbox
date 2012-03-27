function [fx,dF_dX,dF_dTheta] = f_dcm4fmri(Xt,Theta,ut,inF)
% nonlinear DCM evolution function
% function [fx,dF_dX,dF_dTheta] = f_dcm4fmri(Xt,Theta,ut,inF)
% This function evaluates the evolution function of the neuronal states of
% a nonlinear DCM for fMRI model.
%------------------------------------------------------------
% Copyright (C) 2012 Jean Daunizeau / License GNU GPL v2
%------------------------------------------------------------

deltat = inF.deltat;
n = size(Xt,1);

[dxdt,J,dfdp] = get_fx(Xt,Theta,ut,inF);
fx = Xt + deltat.*dxdt;

dF_dX = eye(n) + deltat.*J';
dF_dTheta = deltat.*dfdp;


function [dxdt,J,dfdp] = get_fx(Xt,Theta,ut,inF)

n = size(Xt,1);
nu = size(ut,1);

dfdp = zeros(size(Theta,1),n);
dfdp(inF.indself,:) = -exp(Theta(inF.indself)).*Xt';
In = eye(n);
xI = kron(Xt',In);

A = inF.A;
if isempty(A)
    A = zeros(n,n);
end
indA = inF.indA;
if ~isempty(indA)
    A(A~=0) = Theta(indA);
    A = A - exp(Theta(inF.indself)).*eye(n);
    dfdp(indA,:) = [xI*inF.dA]';
else
    A = A - exp(Theta(inF.indself)).*eye(n);
end

B = inF.B;
indB = inF.indB;
dxB = zeros(n,n);
for i=1:nu
    if isempty(B{i})
        B{i} = zeros(n,n);
    end
    if ~isempty(indB{i})
        B{i}(B{i}~=0) = Theta(indB{i});
        dfdp(indB{i},:) = ut(i).*[xI*inF.dB{i}]';
        dxB = dxB + ut(i).*B{i};
    end
end

C = inF.C;
if isempty(C)
    C = zeros(n,nu);
end
indC = inF.indC;
if ~isempty(indC)
    C(C~=0) = Theta(indC);
    dfdp(indC,:) = [kron(ut',In)*inF.dC]';
end

D = inF.D;
indD = inF.indD;
dxD = zeros(n,n);
dxD2 = dxD;
for i=1:n
    if isempty(D{i})
        D{i} = zeros(n,n);
    end
    if ~isempty(indD{i})
        D{i}(D{i}~=0) = Theta(indD{i});
        tmp = Xt(i)*D{i};
        dxD = dxD + tmp;
        tmp(:,i) = tmp(:,i)+D{i}*Xt;
        dxD2 = dxD2 +tmp;
        dfdp(indD{i},:) = Xt(i)*[xI*inF.dD{i}]';
    end
end

flow = A + dxB + dxD; % Unperturbed flow
dxdt = flow*Xt + C*ut; % vector field
J = A + dxB + dxD2; % Jacobian





