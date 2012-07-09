function [fb] = h_reward_2Q(yt,t,in)
% compares the entry yt with a stored reference answer u0(t)
fb = in.u0(yt+1,t); % yt=0 or 1 (option chosen
