function [gx] = g_wrap_perseveration(x,P,u,in)
% wraps observation function to include (anti-)perseverative bias
% function g_wrapInGame(x,P,u,in)
% A perseveration bias occurs when an agent tends to repeat his past
% action, irrespective of his beliefs and/or preferences. It can be
% modelled as a bonus for repeating the previous move.
% Let V1_t (resp. V0_t) the value of option y_t=1 (resp. y_t=0), at trial
% t. Let R_t be such that: B_t=1 if y_{t-1}=1 and B_t = -1 otherwise. The
% probability P(y_t=1) of picking option y_t=1 is thus given by:
% P(y_t=1) = sig( V1-V0 + delta*B_t)
% where delta is the perseveration weight. Note that if beta<0, then the
% agent tends to alternate his moves.
% IN:
%   - x: hidden states
%   - P: observation params & perseveration weight
%   - u: here, u(2) = agent's last move!
%   - in: includes
%   in.g0 = OBS function handle
%   in.indP0 = indices of params of native OBS function
%   in.indbeta = index of perseveration weight
%   in.in0 = optional input structure to native OBS function
% OUT:
%   - gx: P(y_t=1)

P0 = P(in.indP0); % OBS params of the native OBS function
g0 = feval(in.g0,x,P0,u,in.in0); % P(y_t=1) without persevration

lastMove = u(2);
if VBA_isWeird (lastMove)
    gx = g0; % no perseveration
else
    dV = VBA_sigmoid(g0, 'inverse', true); % V1-V0
    beta = P(in.indbeta); % perseveration weight
    Bt = 2*lastMove - 1; % 1 if y_{t-1}=1 and -1 otherwise 
    gx = VBA_sigmoid(dV + beta.*Bt); % perseverative tendency
end
