function [El] = VBA_ElogBeta(a,b)
% expectation of a log-transformed beta variable
% function [El] = VBA_ElogBeta(a,b)
% IN:
%   - a,b: shape parameters of the beta distribution
% OUT:
%   - El: E[log(x)|a,b]
% Let X be a random variable following a beta distribution with shape
% parameters a and b: X ~ Beta(a,b). This functions derives a
% semi-analytical approximation of  E[log(x)|a,b]. See the technical note
% on this webpage: http://sites.google.com/site/jeandaunizeauswebsite/links/resources.
% NB: the function is meant to be used in cases where b is smaller than 1. 

na = numel(a);
nb = numel(b);
if ~isequal(na,nb)
    disp('Error: vectors of beta shape parameters should have the same size!')
    El = [];
    return
else
    a = VBA_vec(a);
    b = VBA_vec(b);
end

%if any(b<1)
%    warning('beta shape parameter smaller than 1!')
%end
El = zeros(na,1);
ind = find(a>1e-1);
El(ind) = VBA_psi(a(ind)) - VBA_psi(a(ind)+b(ind));
ind = setdiff(1:na,ind);
no = length(ind);
if ~isempty(ind)
    Pf = [
        -1.1659    3.6447  -35.9705
        0.0755   -0.3278   -0.0669
        ];
    P = zeros(no,3);
    for i=1:3
        P(:,i) = g_lu(Pf(:,i),log(b(ind)));
    end
    El(ind) = g_su(P,log(a(ind)));
end

function gx = g_lu(P,u)
gx = P(1) + P(2)*u;

function gx = g_su(P,u)
gx = P(:,3)./(1+exp(-P(:,1).*u+P(:,2)));


