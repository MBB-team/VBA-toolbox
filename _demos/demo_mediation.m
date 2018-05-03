% this demonstrates Baron & Kenny's simple mediation analysis

close all
clear all

n = 32; % sample size
s = 2; % noise standard deviation

X = randn(n,1) + s*10*randn;
M = X + s*randn(n,1) + s*randn;
Y = M + s*randn(n,1) + s*randn;


[out] = mediationAnalysis0(Y,X,M);