function fb = h_goNogo(yt,t,in)
% generates the feedback (go/nogo task) given the agent's last choice
if yt==1 % 'go'
    fb = in.u(t);
else % 'nogo'
    fb = 0;
end
