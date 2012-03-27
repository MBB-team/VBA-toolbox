function  [ fx,dfdx,dfdP ] = f_Qlearn( x_t,P,u_t,in )
%%% Fonction d'apprentissage
% INPUT 
% - x_t : Les états cachés sont les Qvalues (2*1)
% - P : learning rate (1*1)
% - u_t : l'action et la récompense ( 2*1 = 1:actions, 2:récompenses) 
% - in : []
%------------------------------------------------------------
% Copyright (C) 2012 Jean Daunizeau / License GNU GPL v2
%------------------------------------------------------------

alpha = 1./(1+exp(-P));
a = u_t(1)+1;
r = u_t(2);
fx = x_t;
fx(a) = x_t(a) + alpha*(r-x_t(a));

if a == 1
%dfdx = [df1dx1 , df1dx2;
%        df2dx1 , df2dx2]
    dfdx = [1-alpha, 0;
            0, 1];
    dfdP = [alpha*(1-alpha)*(r-x_t(a)),0];
    

elseif a == 2                   
    dfdx = [1, 0;
            0, 1-alpha];
    dfdP = [0,alpha*(1-alpha)*(r-x_t(a))];
    
end

% dfdx = dfdx';
% dfdP = dfdP';