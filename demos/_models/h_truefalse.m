function [fb] = h_truefalse(yt,t,in)
% compares the entry yt with a stored reference answer u0(t)
if isequal(yt,in.u0(t))
    fb = 1;
else
    fb = 0;
end