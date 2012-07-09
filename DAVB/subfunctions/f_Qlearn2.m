function  [ fx,dfdx,dfdP ] = f_Qlearn2( x,P,u,in )
% IN:
% - x_t : Q-values (2*1)
% - P : learning rate (1*1)
% - u_t : previous action and feedback
% - in : []

alpha = 1./(1+exp(-P));
a = u(1)+1; % index of Qvalue to update
r = u(2);
fx = x;
fx(a) = x(a) + alpha*(r-x(a));

if a == 1
    %dfdx = [df1dx1 , df1dx2;
    %        df2dx1 , df2dx2]
    dfdx = [1-alpha, 0;
        0, 1];
    dfdP = [alpha*(1-alpha)*(r-x(a)),0];
    
    
elseif a == 2
    dfdx = [1, 0;
        0, 1-alpha];
    dfdP = [0,alpha*(1-alpha)*(r-x(a))];
    
end
    
% dfdx = dfdx';
% dfdP = dfdP';