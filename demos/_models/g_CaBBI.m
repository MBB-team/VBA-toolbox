% The observation function:
%               
%               function [F] = g_CaBBI(Xt,Phi,I_inp,inG)
%
% this function maps [Ca2+] kinetics to fluorescence observations (F) 
% through a non-linear saturating function


function [F] = g_CaBBI(Xt,Phi,I_inp,inG)

ind = inG.ind;                 
Ca = Xt(ind);                        % [Ca2+] kinetics
d_F =  Phi(1);                       % offset parameter
Kd =  200;                           % dissociation constant
scale = inG.k_F0*exp(Phi(2));        % scale paraemeter
F =  scale*(Ca/(Ca + Kd)) + d_F;     % nonlinear mapping from [Ca2+] to F

end
