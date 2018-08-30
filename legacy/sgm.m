function y = sgm(x,a)
% legacy code
s = warning ('on');
warning ('*** The function `sgm` is now deprecated. Please see `VBA_sigmoid` for an alternative.') 
warning (s);

% fallback
if nargin == 1
    a = 1;
end
y = VBA_sigmoid (x, 'scale', a);
