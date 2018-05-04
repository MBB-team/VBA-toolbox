function [fx,dfdx,dfdP] = logExp(x,P)
% approximates the mappings y=a at x=-Inf and y=b(x-c) at x=Inf 
% [fx,dfdx,dfdP] = logExp(x,P)
% IN:
%   - x: the variable to be passed through the logExp mapping
%   - P: parameters of the logExp mapping (P = [a b c])
% OUT:
%   - fx: the logExp mapping
%   - dfdx: the gradient wrt to x
%   - dfdP: the gradient wrt to P

try;a=P(1);catch;a=0;end % offset at x=-Inf
try;b=P(2);catch;b=1;end % linear slope at x=Inf
try;c=P(3);catch;c=0;end % elbow position

x = x(:)';
ea = exp(a);
ebx = exp(b.*(x-c));
fx = log(ea+ebx);
sx = ea./(ea+ebx);
dfdx = b.*(1 - sx);
dfdP = [ sx ; (x-c).*(1-sx) ; -b.*(1 - sx) ];