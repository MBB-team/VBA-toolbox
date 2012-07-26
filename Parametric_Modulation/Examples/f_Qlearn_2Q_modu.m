function  [ fx] = f_Qlearn_2Q_modu( x_t,P,u_t,in )
%,dfdx,dfdP
%%% update of Qvalue a scalar value given past Qvalue and reward
% INPUT 
% - x_t : Qvalue at trial t-1
% - P : learning rate
% - u_t : 
%    - Action performed (binary)
%    - Obtained Reward 


% Transformation on parameters
P(1) = P(1) + P(2)*u_t(3) + P(3)*u_t(4) ; % modulation of alpha


try 
    if isequal(in.param_transform.type,'sigmoid')
        alpha = 1./(1+exp(-P(1)));
    elseif isequal(in.param_transform.type,'modified sigmoid')
        alpha = in.param_transform.a + (in.param_transform.b-in.param_transform.a)./(1+exp(-P(1)));                
    end
catch
    alpha = P(1);
end


a = u_t(1)+1; % index of Qvalue to update (1 or 2)
r = u_t(2); % 2!
fx = x_t;
fx(a) = x_t(a) + alpha*(r-x_t(a));


