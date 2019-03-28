function [fx,dF_dX,dF_dTheta] = f_dcm_extension(Xt,Theta,ut,inF)
% nonlinear DCM evolution function
% function [fx,dF_dX,dF_dTheta] = f_dcm_extension(Xt,Theta,ut,inF)
% This function evaluates the evolution function of the neuronal states of
% a nonlinear DCM for fMRI model.

deltat = inF.deltat;
n = size(Xt,1);
nr = length(inF.r);

[dxdt,J,dfdp] = get_fx(Xt,Theta,ut,inF);

fx = Xt(end-nr+1:end) + deltat.*dxdt;
dF_dX = [zeros(n-nr,nr) ; eye(nr)] + deltat.*J';
dF_dTheta = deltat.*dfdp;


function [dxdt,J,dfdp] = get_fx(Xt,Theta,ut,inF)


nr = length(inF.r) ;
n = size(Xt,1) - nr;
%n=nz+nr;
nu = size(ut,1);
idR =(1:nr) + n;

% == self
hself = - diag(exp(Theta(inF.indhself))) ;

xI = easykron(Xt(1:n),n,nr);

%%

dfdp = zeros(size(Theta,1),nr);
dfdp(inF.indhself,:) = diag(diag(hself).*Xt(idR));
%dfdp(inF.indconst,:) = diag(ones(1,nr));

%%

%== A
hA = inF.hA;

indhA = inF.indhA;
if ~isempty(indhA)
    hA(hA~=0) = Theta(indhA);
    dfdp(indhA,:) = inF.dhA'*xI;
end

%%== B
hB = inF.hB;
indhB = inF.indhB;
dxhB = zeros(nr,n);
for i=1:nu
    if ~isempty(indhB{i})
        hB{i}(hB{i}~=0) = Theta(indhB{i});
        dfdp(indhB{i},:) = ut(i).*(inF.dhB{i}'*xI);
        dxhB = dxhB + ut(i).*hB{i};
    end
end

%== C
hC = inF.hC;
indhC = inF.indhC;
if ~isempty(indhC)
    hC(hC~=0) = Theta(indhC);
    dfdp(indhC,:) = inF.dhC'*easykron(ut,nu,nr);
end

%== D
hD = inF.hD;
indhD = inF.indhD;
dxhD = zeros(nr,n);
dxhD2 = dxhD;
for i=1:n
    if ~isempty(indhD{i})
        hD{i}(hD{i}~=0) = Theta(indhD{i});
        tmp = Xt(i)*hD{i};
        dxhD = dxhD + tmp;
        tmp(:,i) = tmp(:,i)+hD{i}*Xt(1:n);
        dxhD2 = dxhD2 +tmp;
        dfdp(indhD{i},:) = Xt(i)*(inF.dhD{i}'*xI);
    end
end

%consts = Theta(inF.indconst);

%%


flow = hA + dxhB + dxhD; % Unperturbed flow
dxdt = hself*Xt(idR) + flow*Xt(1:n) + hC*ut ; %+ consts; % vector field
Jself = [zeros(n,nr) ; hself]';
J = [hA + dxhB + dxhD2 , zeros(nr,nr)] + Jself ; % Jacobian

function B=easykron(X,nx,n)
B = zeros(n*nx,n);
  idx = 1:n:(n*nx+eps);
  for k = 1:n  
    B(idx,k) = X;
    idx = idx+1;
  end




