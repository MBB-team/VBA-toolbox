function [fx,dF_dX,dF_dTheta] = f_dcm4fmri(Xt,Theta,ut,inF)
% nonlinear DCM evolution function
% function [fx,dF_dX,dF_dTheta] = f_dcm4fmri(Xt,Theta,ut,inF)
% This function evaluates the evolution function of the neuronal states of
% a nonlinear DCM for fMRI model.

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

xI = easykron(Xt,n,n);

A = inF.A;
indA = inF.indA;
if ~isempty(indA)
    A(A~=0) = Theta(indA);
    A = A - exp(Theta(inF.indself)).*eye(n);
    dfdp(indA,:) = inF.dA'*xI;
else
    A = A - exp(Theta(inF.indself)).*eye(n);
end

B = inF.B;
indB = inF.indB;
dxB = zeros(n,n);
for i=1:nu
%     if isempty(B{i})
%         B{i} = zeros(n,n);
%     end
    if ~isempty(indB{i})
        B{i}(B{i}~=0) = Theta(indB{i});
        dfdp(indB{i},:) = ut(i).*(inF.dB{i}'*xI);
        dxB = dxB + ut(i).*B{i};
    end
end

C = inF.C;
% if isempty(C)
%     C = zeros(n,nu);
% end
indC = inF.indC;
if ~isempty(indC)
    C(C~=0) = Theta(indC);
    dfdp(indC,:) = inF.dC'*easykron(ut,nu,n);
end

D = inF.D;
indD = inF.indD;
dxD = zeros(n,n);
dxD2 = dxD;
for i=1:n
    if ~isempty(indD{i})
        D{i}(D{i}~=0) = Theta(indD{i});
        tmp = Xt(i)*D{i};
        dxD = dxD + tmp;
        tmp(:,i) = tmp(:,i)+D{i}*Xt;
        dxD2 = dxD2 +tmp;
        dfdp(indD{i},:) = Xt(i)*(inF.dD{i}'*xI);
    end
end

flow = A + dxB + dxD; % Unperturbed flow
dxdt = flow*Xt + C*ut; % vector field
J = A + dxB + dxD2; % Jacobian

function B=easykron(X,nx,n)
  B = zeros(n*nx,n);
  idx = 1:n:(n*nx+eps);
  for k = 1:n 
    B(idx,k) = X;
    idx = idx+1;
  end





