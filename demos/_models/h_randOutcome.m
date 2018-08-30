function [fb] = h_randOutcome(yt,t,in)
% compares the entry yt with a stored reference answer u0(t)
if isequal(yt,in.u0(t))
    fb = in.er;
else
    fb = -in.er;
end
fb = fb + sqrt(in.vr)*randn;