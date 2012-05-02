function  [ fx,dfdx,dfdP ] = f_Qlearn( x_t,P,u_t,in )
% IN: 
% - x_t : Q-values (2*1)
% - P : learning rate (1*1)
% - u_t : previous action and feedback
% - in : []

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