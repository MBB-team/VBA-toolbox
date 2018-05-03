function [g,dgdx,dgdP] = g_logistic(x,P,u,in)
% derives the probability of outcome variable y=1, under the logistic model
g = sig(u'*P);
dgdx = [];
dgdP = u*diag(g.*(1-g));



function s = sig(x)
s = 1./(1+exp(-x));
% s(s<1e-4) = 1e-4;
% s(s>1-1e-4) = 1-1e-4;

