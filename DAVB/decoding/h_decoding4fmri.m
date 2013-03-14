function [hx,dH_dR,dH_dPsi] = h_decoding4fmri(Rt,Psi,XUt,inH)
% nonlinear DCM evolution function
% function [hx,dH_dX,dH_dPsi] = h_decoding4fmri(Rt,Psi,XUt,inH)
% This function evaluates the nonlinear deconding function of the neuronal states of
% DCM for fMRI model.
% - IN
%      Xt:  actual state of the system
%      Psi: parameters of the decoding scheme
%      ut:  actual inputs
%      inH: summary structure of the decoding scheme
% - OUT
%      hx: state at the next timestep
%      dH_dX,dH_dPsi: instantaneus partial derivative of the function with respect to
%                     state (dynamic response prediction) and parameter respectively

deltat = inH.deltat;

% compute derivatives
[drdt,J,dhdp] = get_hx(Rt,Psi,XUt,inH);

% euler integration
hx = Rt + deltat.*drdt;

% partial derivatives
dH_dR = eye(inH.n_r) + deltat.*J';
dH_dPsi = deltat.*dhdp;

end

function [drdt,J,dhdp] = get_hx(Rt,Psi,XUt,inH)

% size of the state, input and response vectors
Xt = XUt(inH.nx) ;
ut = XUt(inH.nu) ;
n = size(Xt,1);
nu = size(ut,1);
nr = inH.n_r;

% == recover matrices of the taylor decomposition
% direct input prediction
hU = inH.hU;
indhU = inH.indhU ;
if isempty(hU)
    hU = zeros(nr,nu);
end
if ~isempty(indhU)
    indhU = inH.indhU;
    hU(hU~=0) = Psi(indhU);  
end
% input interaction
hU2 = inH.hU2;
indhU2 = inH.indhU2;
for i=1:nu
    if isempty(hU2{i})
        hU2{i} = zeros(nr,nu);
    end
    if ~isempty(indhU2{i})
        hU2{i}(hU2{i}~=0) = Psi(indhU2{i});
    end
end
% direct state prediction
hX = inH.hX;
if isempty(hX)
    hX = zeros(nr,n);
else
    indhX = inH.indhX;
    hX(hX~=0) = Psi(indhX);  
end
% state interaction
hX2 = inH.hX2;
indhX2 = inH.indhX2;
for i=1:n
    if isempty(hX2{i})
        hX2{i} = zeros(nr,n);
    end
    if ~isempty(indhX2{i})
        hX2{i}(hX2{i}~=0) = Psi(indhX2{i});
    end
end
% bilinear interactions
hBi = inH.hBi;
indhBi = inH.indhBi;
for i=1:nu
    if isempty(hBi{i})
        hBi{i} = zeros(nr,n);
    end
    if ~isempty(indhU2{i})
        hBi{i}(hBi{i}~=0) = Psi(indhBi{i});
    end
end
% self deactivation
hself = - exp(Psi(inH.indhself)).*eye(nr) ;

% == computing derivatives

% time derivative
dhX2 = zeros(nr,n);
dhX2J = zeros(nr,n);
for i=1:n 
    tmp = Xt(i)*hX2{i};
    dhX2 = dhX2 + tmp ;
    tmp(:,i) = tmp(:,i)+hX2{i}*Xt;
    dhX2J = dhX2J +tmp;      
end
dhBi = zeros(nr,n);
dhU2 = zeros(nr,nu);
for i=1:nu
    dhBi = dhBi + ut(i).*hBi{i};
    dhU2 = dhU2 + ut(i).*hU2{i};
end

decoding = hX + dhX2 + dhBi ;
computational = hU + dhU2 ;
drdt = decoding*Xt + computational*ut + hself*Rt;


% jacobian
J = hself ; 

% derivative w.r. parameters
In = eye(nr);
xI = kron(Xt',In);
uI = kron(ut',In);
dhdp = zeros(numel(Psi),nr);
% hX
if ~isempty(inH.indhX)
    dhdp(indhX,:) = [xI*inH.dhX]'; 
end
% hX2
for i=1:n
    if ~isempty(inH.indhX2{i})
        dhdp(indhX2{i},:) = Xt(i)*[xI*inH.dhX2{i}]';
    end
end
% hU
if ~isempty(inH.indhU)
    dhdp(indhU,:) = [uI*inH.dhU]'; 
end
% hU2, hBi
for i=1:nu
    if ~isempty(inH.indhU2{i})
     dhdp(indhU2{i},:) = ut(i)*[uI*inH.dhU2{i}]';
    end
    if  ~isempty(inH.indhBi{i})
        dhdp(indhBi{i},:) = ut(i).*[xI*inH.dhBi{i}]';
    end
end
% self
dhdp(inH.indhself,:) = -exp(Psi(inH.indhself)).*Rt';

end








