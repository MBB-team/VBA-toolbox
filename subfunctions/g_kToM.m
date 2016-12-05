function [gx] = g_kToM(x,P,u,in)
% wrapper around recursive k-ToM observation function
gx = ObsRecGen(x,P,u,in);