function  [fx] = f_Qlearn_gammaLR(x,P,u,in)
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

gamma = exp(P(1));
tau = exp(P(2));
alpha = 1./(1+exp(-x(3))); % learning rate is bounded between 0 and 1.
fx = zeros(4,1);
fx(1) = x(1) + alpha*(u(2)-x(1))*u(1);
fx(2) = x(2) + alpha*(u(2)-x(2))*(1-u(1));
fx(3) = x(3) + x(4);
fx(4) = x(4) + gamma.*tau.*u(3) - 2*tau.*x(4) - tau^2.*x(3);
