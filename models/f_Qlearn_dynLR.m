function  [fx,dfdx] = f_Qlearn_dynLR(x,P,u,in)
% evolution function of q-values of a RL agent with dynamic learning rate
% [fx] = f_Qlearn_dynLR(x,P,u,in)
% Here, there are only two q-values to evolve, i.e. there are only two
% actions to reinforce (2-armed bandit problem).
% IN:
%   - x_t : q-values and learning rate
%   - P : [useless]
%   - u : u(1)=previous action, u(2)=feedback
%   - in : [useless]
% OUT:
%   - fx: evolved q-values and learning rate (3x1)

alpha = 1./(1+exp(-x(3))); % learning rate is bounded between 0 and 1.
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
    dfdx = [1-alpha, 0, 0;
        0, 1, 0;
        alpha*(1-alpha)*(r-x(a)), 0, 1];
elseif a == 2
    dfdx = [1, 0, 0;
        0, 1-alpha, 0;
        0, alpha*(1-alpha)*(r-x(a)), 1];
end