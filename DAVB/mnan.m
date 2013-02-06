function y = mnan(x,dim)
% mean operator, having removed NaN entries
xin = ones(size(x));
xin(isnan(x)) = 0;
if nargin==1, 
  dim = find(size(x)~=1, 1 );
  if isempty(dim), dim = 1; end
  y = snan(x)./sum(xin,dim);
else
  y = snan(x,dim)./sum(xin,dim);
end