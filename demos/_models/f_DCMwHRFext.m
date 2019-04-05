function [fx,dfdx,dfdp] = f_DCMwHRFext(Xt,Theta,ut,inF) % 
% DCM for fMRI evolution function (including Balloon model)
% function [fx,dF_dX,dF_dTheta] = f_DCMwHRF(Xt,Theta,ut,inF)
% This function evaluates the evolution function DCM for fMRI, including
% the Balloon HRF model.


ut = ut(inF.confounds.indu);
nx = length(Xt);
nr = length(inF.r);

%- hidden states evolution
xn = Xt(inF.n5);
if ~isfield(inF,'fast')
    [fxh,dfdxh,dfdph] = f_HRF3(Xt(1:nx-nr),Theta,xn,inF);
else
    fxh = zeros(nx-nr,1);
    dfdxh = zeros(nx-nr,nx-nr);
    dfdph = zeros(length(Theta),nx-nr);
end

[fxn,dfdxn,dfdpn] = f_dcm4fmri(Xt(inF.n5),Theta,ut,inF);
[fxr,dfdxr,dfdpr] = f_dcm_extension(Xt([inF.n5 inF.r]),Theta,ut,inF);

%== Reshape flow field and gradients
%- flow
fx = zeros(nx,1);
fx(1:nx-nr) = fxh ;
fx(inF.n5) = fxn ;
fx(inF.r) = fxr ;

%- jacobian
dfdx = zeros(nx,nx);
dfdx(1:nx-nr,1:nx-nr) = dfdxh ;
dfdx(inF.n5,inF.n5) = dfdxn;
dfdx([inF.n5 inF.r],inF.r) = dfdxr ; 

%- wrt parameters
dfdp = zeros(length(Theta),nx);
dfdp(:,1:nx-nr) = dfdph ;
dfdp(:,inF.n5) = dfdp(:,inF.n5)+dfdpn ;
dfdp(:,inF.r) = dfdpr ;



end





