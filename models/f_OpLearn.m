function  fx = f_OpLearn( x,P,u,in )
% re-directs VB updates to the previously chosen action
% IN:
% - x : 1st and 2nd order moments of the action-dependent beliefs
% - P : observer's prior parameters
% - u_t : previous action and feedback
% - in : []

a = u(1);%+1; % action that was chosen on the previous trial
r = u(2); % feedback received on the previous trial

if a==1
    % update beliefs conditional on first action
    ixu = 1:5;
else
    % update beliefs conditional on second action
    ixu = 6:10;
end

fx = x;
fx(ixu) = f_VBvolatile0(x(ixu),P,r,in); % learning from outcome category (0 or 1), not from reward value


