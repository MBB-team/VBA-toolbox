function [gx] = g_wk2004(x,P,u,in)
% Function to perform optiomal prior and cue combination to correct behavior
% as described in 'Bayesian integration in sensorimotor learning'

% INPUT
% - x : [] empty
% - P : 
%   - [1] mean of prior on deviation 
%   - [2] std  of prior on deviation
% - u : [] empty
% - in : 
%   - [1] mean of mid course flash (also real shift)
%   - [2] std  of mid course flash
% OUTPUT
% - gx : actual correction


% initial position : 0
x_s = in.shifts;
s_s = in.std_cues;
mu_p = P(1);
s_p  = exp(P(2));

gx = s_s.^2./(s_s.^2+s_p^2)*mu_p + s_p^2./(s_s.^2+s_p^2).*x_s; % predicted shift
gx=gx';
end