function   y = sigmoid(x)
% legacy code
s = warning ('on');
warning ('*** The function `sigmoid` is now deprecated. Please see `VBA_sigmoid` for an alternative.') 
warning (s);

% fallback
y = VBA_sigmoid (x);