function ux = u_prospect_theory(x,P,in)
% computes the utility of an outcome x'.
% IN:
%   - x: outcome.
%   - P: parameters of the utility function.
%   - in: further quantities handed to the function.
% OUT:
%   - ux : utility of the outcome

Ip = find(x>=0); Im = find(x<0);
ux = zeros(size(x));
ux(Ip) = exp(P(1))*log(x(Ip)+1);
ux(Im) = -exp(-P(1))*log(-x(Im)+1);
