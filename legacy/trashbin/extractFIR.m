function out = extractFIR(x,n)
out = reshape(VBA_vec(x)',n,[]);
out = mean(out,2)';
