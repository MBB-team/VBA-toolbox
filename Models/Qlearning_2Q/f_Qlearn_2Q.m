function  [ fx] = f_Qlearn_2Q( x_t,P,u_t,in )
%,dfdx,dfdP
%%% update of Qvalue a scalar value given past Qvalue and reward
% INPUT 
% - x_t : Qvalue at trial t-1
% - P : learning rate
% - u_t : 
%    - Action performed (binary)
%    - Obtained Reward 


% Transformation on parameters
try 
    if isequal(in.param_transform.type,'sigmoid')
        alpha = 1./(1+exp(-P));
    elseif isequal(in.param_transform.type,'modified sigmoid')
        alpha = in.param_transform.a + (in.param_transform.b-in.param_transform.a)./(1+exp(-P));                
    end
catch
    alpha = P;
end

a = u_t(1)+1; % index of Qvalue to update (1 or 2)
r = u_t(2); % 2!
fx = x_t;
fx(a) = x_t(a) + alpha*(r-x_t(a));

if a == 1
    dfdx = [1-alpha, 0;
        0, 1];
    dfdP = [alpha*(1-alpha)*(r-x_t(a)),0];
    
    
elseif a == 2
    dfdx = [1, 0;
        0, 1-alpha];
    dfdP = [0,alpha*(1-alpha)*(r-x_t(a))];
    
end

 %fx = fx*nan; (to test wrong outputs!)
