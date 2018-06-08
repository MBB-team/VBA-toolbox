function y = num2multi(n,s1)

if ~exist('s1')
    s1 = max(n);
end

n = VBA_vec(n)';
s2 = size(n,2);

y = +(repmat(n,s1,1) == repmat((1:s1)',1,s2)) ;

