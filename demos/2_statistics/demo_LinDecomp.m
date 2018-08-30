% demo for linear decomposition

clear all
close all

k = 5; % number of simulated components
n = 15; % 1st dim of data matrix
p = 7; % 2nd dim of data matrix

A0 = randn(n,k);
B0 = randn(k,p);
Y = A0*B0;
y = Y+randn(size(Y));

[A,B,out,posterior] = VBA_LinDecomp(y,k);


hf = figure('color',[1 1 1]);
ha = subplot(2,2,1,'parent',hf);
plot(ha,A,A0,'.')
ha = subplot(2,2,2,'parent',hf);
plot(ha,B',B0','.')
ha = subplot(2,2,3,'parent',hf);
bar(corr(A,A0),'parent',ha)
ha = subplot(2,2,4,'parent',hf);
bar(corr(B',B0'),'parent',ha)