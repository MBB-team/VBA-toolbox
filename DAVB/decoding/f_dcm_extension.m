function [fx,dF_dX,dF_dTheta] = f_dcm_extension(Xt,Theta,ut,inF)
% nonlinear DCM evolution function
% function [fx,dF_dX,dF_dTheta] = f_dcm4fmri_ext(Xt,Theta,ut,inF)
% This function evaluates the evolution function of the neuronal states of
% a nonlinear DCM for fMRI model.

deltat = inF.deltat;
n = size(Xt,1);
nr = length(inF.r);

[dxdt,J,dfdp] = get_fx(Xt,Theta,ut,inF);
fx = Xt(end-nr+1:end) + deltat.*dxdt;


dF_dX = [zeros(n-nr,nr) ; eye(nr)] + deltat.*J';
% dF_dX=dF_dX'; % TODO: transpose transpose etc.

dF_dTheta = deltat.*dfdp;


function [dxdt,J,dfdp] = get_fx(Xt,Theta,ut,inF)

nr = length(inF.r) ;
n = size(Xt,1) - nr;
%n=nz+nr;
nu = size(ut,1);
idR =(1:nr) + n;

dfdp = zeros(size(Theta,1),nr);
dfdp(inF.indhself,:) = -exp(Theta(inF.indhself)).*Xt(idR)';
dfdp(inF.indconst,:) = diag(ones(1,nr));

In = eye(nr);
xI = kron(Xt(1:n)',In);

% == self
hself = - exp(Theta(inF.indhself)).*eye(nr);

%== A
hA = inF.hA;
if isempty(hA)
    hA = zeros(nr,n);
end
indhA = inF.indhA;
if ~isempty(indhA)
    hA(hA~=0) = Theta(indhA);
    dfdp(indhA,:) = [xI*inF.dhA]';
end

%%== B
hB = inF.hB;
indhB = inF.indhB;
dxhB = zeros(nr,n);
for i=1:nu
    if isempty(hB{i})
        hB{i} = zeros(nr,n);
    end
    if ~isempty(indhB{i})
        hB{i}(hB{i}~=0) = Theta(indhB{i});
        dfdp(indhB{i},:) = ut(i).*[xI*inF.dhB{i}]';
        dxhB = dxhB + ut(i).*hB{i};
    end
end

%== C
hC = inF.hC;
if isempty(hC)
    hC = zeros(nr,nu);
end
indhC = inF.indhC;
if ~isempty(indhC)
    hC(hC~=0) = Theta(indhC);
    dfdp(indhC,:) = [kron(ut',In)*inF.dhC]';
end

%== D
hD = inF.hD;
indhD = inF.indhD;
dxhD = zeros(nr,n);
dxhD2 = dxhD;
for i=1:n
    if isempty(hD{i})
        hD{i} = zeros(nr,n);
    end
    if ~isempty(indhD{i})
        hD{i}(hD{i}~=0) = Theta(indhD{i});
        tmp = Xt(i)*hD{i};
        dxhD = dxhD + tmp;
        tmp(:,i) = tmp(:,i)+hD{i}*Xt(1:n);
        dxhD2 = dxhD2 +tmp;
        dfdp(indhD{i},:) = Xt(i)*[xI*inF.dhD{i}]';
    end
end

% consts = zeros(n,1);
consts = Theta(inF.indconst);

%%


flow = hA + dxhB + dxhD; % Unperturbed flow
dxdt = hself*Xt(idR) + flow*Xt(1:n) + hC*ut + consts; % vector field
Jself = [zeros(n,nr) ; hself]';
J = [hA + dxhB + dxhD2 , zeros(nr,nr)] + Jself ; % Jacobian





