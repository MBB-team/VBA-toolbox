function [fx] = f_kToM(x,P,u,inF)
% wrapper around recursive k-ToM evolution function
[fx,indlev] =RecToMfunction(x,P,u,inF);