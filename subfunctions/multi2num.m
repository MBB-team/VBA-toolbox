function n = multi2num(y)


[s1 s2] = size(y);

n = max(y .* repmat((1:s1)',1,s2)) ;

