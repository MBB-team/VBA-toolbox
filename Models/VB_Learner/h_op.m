function [fb] = h_op(yt,t,in)
% in.u0 is 2 dimensional :
% it contains the two possible outputs of the decision
fb = in.u0(yt+1,t);
end