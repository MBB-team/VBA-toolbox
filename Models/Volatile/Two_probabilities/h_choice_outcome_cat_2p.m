function [fb] = h_choice_outcome_cat_2p(yt,t,in)
% returns the value of the outcome cat of a binary choice
% IN
% - yt : the binary decision (0 or 1)
% OUT
% - fb : the binary outcome for decision yt (0 or 1) 
fb = in.u0(yt+1,t);
