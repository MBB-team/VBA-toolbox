function out = extractFIR(x,n)
out = reshape(vec(x)',n,[]);
out = mean(out,2)';
