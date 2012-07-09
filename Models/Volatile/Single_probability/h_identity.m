function [fb] = h_identity(yt,t,in)
% feedback is entry (no really need for feedback here...)
fb = in.u0(t);
