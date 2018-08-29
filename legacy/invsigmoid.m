function  y = invsigmoid(x)
% legacy code
s = warning ('on');
warning ('*** The function `invsigmoid` is now deprecated. Please see `VBA_sigmoid` for an alternative.') 
warning (s);

% fallback
y = VBA_sigmoid (x, 'inverse', true);