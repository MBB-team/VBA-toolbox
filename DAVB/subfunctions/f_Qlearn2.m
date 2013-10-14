function  [fx,dfdx,dfdP] = f_Qlearn2(x,P,u,in)
% evolution function of q-values of a RL agent (2-armed bandit problem)
% [fx,dfdx,dfdP] = f_Qlearn2(x,P,u,in)
% Here, there are only two q-values to evolve, i.e. there are only two
% actions to reinforce (2-armed bandit problem).
% IN:
%   - x_t : q-values (2x1)
%   - P : (inverse-sigmoid) learning-rate
%   - u : u(1)=previous action, u(2)=feedback
%   - in : [useless]
% OUT:
%   - fx: evolved q-values (2x1)
%   - dfdx/dfdP: gradient of the q-values evolution function, wrt q-avlues
%   and evolution parameter, respectively.

alpha = 1./(1+exp(-P)); % learning rate is bounded between 0 and 1.
r = u(2);
fx = x;

% index of Qvalue to update
if u(1)==1
    a = 1;
else
    a = 2;
end
fx(a) = x(a) + alpha*(r-x(a));

if a == 1
    dfdx = [1-alpha, 0;
        0, 1];
    dfdP = [alpha*(1-alpha)*(r-x(a)),0];
elseif a == 2
    dfdx = [1, 0;
        0, 1-alpha];
    dfdP = [0,alpha*(1-alpha)*(r-x(a))];
end