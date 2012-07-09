function  fx = f_OpLearn_2p( x,P,u,in )
% IN:
% - x : 1st and 2nd order moments of the action-dependent beliefs
% - P : observer's prior parameters
% - u_t : previous 
%       - selected option (binary)
%       - obtained outcome category (binary)
% - in : []

a = u(1)+1; % action that was selected on the previous trial
r = u(2); % binary outcome category for selected action

if a==1
    % update beliefs conditional on first action
    ixu = 1:5;
elseif a==2
    % update beliefs conditional on second action
    ixu = 6:10;
end

fx = x;
fx(ixu) = f_VBvolatile_1p(x(ixu),P,u(2),in); % learning from outcome category (0 or 1), not from reward value


