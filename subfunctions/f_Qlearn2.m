function  [fx,dfdx,dfdP] = f_Qlearn2(x,P,u,in)
% evolution function of q-values of a RL agent (2-armed bandit problem)
% [fx,dfdx,dfdP] = f_Qlearn2(x,P,u,in)
% Here, there are only two q-values to evolve, i.e. there are only two
% actions to reinforce (2-armed bandit problem).
% IN:
%   - x_t : q-values (2x1)
%   - P : (inverse-sigmoid) learning-rate
%   - u : u(1)=previous action (1 or 0), u(2)=feedback
%   - in : [useless]
% OUT:
%   - fx: evolved q-values (2x1)
%   - dfdx/dfdP: gradient of the q-values evolution function, wrt q-avlues
%   and evolution parameter, respectively.

alpha = 1./(1+exp(-P)); % learning rate is bounded between 0 and 1.
fx = zeros(2,1);
pe = u(2)-x; % prediction error
fx(1) = x(1) + alpha*pe(1)*u(1);
fx(2) = x(2) + alpha*pe(2)*(1-u(1));
% gradients' derivation
if u(1)==1
    dfdx = [1-alpha, 0;
            0, 1];
    dfdP = [alpha*(1-alpha)*pe(1),0];
else
    dfdx = [1, 0;
            0, 1-alpha];
    dfdP = [0,alpha*(1-alpha)*pe(2)];
end