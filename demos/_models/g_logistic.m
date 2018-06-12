function [g,dgdx,dgdP] = g_logistic(x,P,u,in)
% derives the probability of outcome variable y=1, under the logistic model
g = VBA_sigmoid(u'*P);
dgdx = [];
dgdP = u*diag(g.*(1-g));