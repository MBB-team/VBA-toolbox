function y = snan(x,dim)
% sum operator, having removed NaN entries
x(isnan(x)) = 0;
if nargin==1, 
  dim = find(size(x)~=1, 1 );
  if isempty(dim), dim = 1; end
  y = sum(x);
else
  y = sum(x,dim);
end