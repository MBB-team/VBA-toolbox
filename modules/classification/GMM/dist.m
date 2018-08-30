function d = dist(y)
n = size(y,2);
d = zeros(n,n);
for i=1:n
    tmp = repmat(y(:,i),1,n) - y;
    tmp = tmp.^2;
    d(i,:) = sum(tmp,1);
end
