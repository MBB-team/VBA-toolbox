function [fx,dF_dX,dF_dTheta] = f_dcm4fmri_ext(Xt,Theta,ut,inF)
% nonlinear DCM evolution function
% function [fx,dF_dX,dF_dTheta] = f_dcm4fmri_ext(Xt,Theta,ut,inF)
% This function evaluates the evolution function of the neuronal states of
% a nonlinear DCM for fMRI model.

deltat = inF.deltat;
n = size(Xt,1);

[dxdt,J,dfdp] = get_fx(Xt,Theta,ut,inF);
fx = Xt + deltat.*dxdt;

dF_dX = eye(n) + deltat.*J';
dF_dTheta = deltat.*dfdp;


function [dxdt,J,dfdp] = get_fx(Xt,Theta,ut,inF)

nr = length(inF.r) ;
nz = size(Xt,1) - nr;
n=nz+nr;
nu = size(ut,1);
idZ = 1:nz;
idR =(1:nr) + nz;

dfdp = zeros(size(Theta,1),n);
dfdp(inF.indself,idZ) = -exp(Theta(inF.indself)).*Xt(idZ)';
dfdp(inF.indhself,idR) = -exp(Theta(inF.indhself)).*Xt(idR)';
dfdp(inF.indconst,idR) = ones(1,nr);

In = eye(n);
xI = kron(Xt',In);

%== A
A = inF.A;
if isempty(A)
    A = zeros(n,n);
end
indA = inF.indA;
if ~isempty(indA)
    A(A~=0) = Theta(indA);
    A(idZ,idZ) = A(idZ,idZ) - exp(Theta(inF.indself)).*eye(nz);
    dfdp(indA,:) = [xI*inF.dA]';
else
    A = A - exp(Theta(inF.indself)).*eye(n);
end

hA = inF.hA;
if isempty(hA)
    hA = zeros(n,n);
end
indhA = inF.indhA;
if ~isempty(indhA)
    hA(hA~=0) = Theta(indhA);
    hA(idR,idR) = hA(idR,idR) - exp(Theta(inF.indhself)).*eye(nr);
    dfdp(indhA,:) = [xI*inF.dhA]';
else
    hA(idR,idR) = hA(idR,idR) - exp(Theta(inF.indhself)).*eye(nr);
end

%%== B
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

hB = inF.hB;
indhB = inF.indhB;
dxhB = zeros(n,n);
for i=1:nu
    if isempty(hB{i})
        hB{i} = zeros(n,n);
    end
    if ~isempty(indhB{i})
        hB{i}(hB{i}~=0) = Theta(indhB{i});
        dfdp(indhB{i},:) = ut(i).*[xI*inF.dhB{i}]';
        dxhB = dxhB + ut(i).*hB{i};
    end
end

%== C
C = inF.C;
if isempty(C)
    C = zeros(n,nu);
end
indC = inF.indC;
if ~isempty(indC)
    C(C~=0) = Theta(indC);
    dfdp(indC,:) = [kron(ut',In)*inF.dC]';
end

hC = inF.hC;
if isempty(C)
    hC = zeros(n,nu);
end
indhC = inF.indhC;
if ~isempty(indhC)
    hC(hC~=0) = Theta(indhC);
    dfdp(indhC,:) = [kron(ut',In)*inF.dhC]';
end

%== D
D = inF.D;
indD = inF.indD;
dxD = zeros(n,n);
dxD2 = dxD;
for i=1:nz
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

hD = inF.hD;
indhD = inF.indhD;
dxhD = zeros(n,n);
dxhD2 = dxhD;
for i=1:nz
    if isempty(hD{i})
        hD{i} = zeros(n,n);
    end
    if ~isempty(indD{i})
        hD{i}(hD{i}~=0) = Theta(indhD{i});
        tmp = Xt(i)*hD{i};
        dxhD = dxhD + tmp;
        tmp(:,i) = tmp(:,i)+hD{i}*Xt;
        dxhD2 = dxhD2 +tmp;
        dfdp(indhD{i},:) = Xt(i)*[xI*inF.dhD{i}]';
    end
end

%%

A = A+hA;
dxB = dxB+dxhB;
C = C + hC ;
dxD = dxD+dxhD;
dxD2 = dxD2+dxhD2;

flow = A + dxB + dxD; % Unperturbed flow
dxdt = flow*Xt + C*ut; % vector field
J = A + dxB + dxD2; % Jacobian





