function [y] = medianfilter0(x,K,verbose)
% operates a simple median filter (on the second dimension of x)
% function [y] = medianfilter0(x,K)
% In brief, the function slides a window of size K, and replaces the signal
% at the middle of the window by the median signal within the window, i.e.:
% y(:,i) = median(x(:,i-K/2:i+K/2),2)
% IN:
%   - x: a pXn matrix
%   - K: the size of the sliding window
%   - verbose: flag for verbose mode {default=1}
% OUT:
%   - y: themedian-filtered signal

try,verbose;catch,verbose=1;end
[p,n] = size(x);
x = [repmat(x(:,1),p,ceil(K/2)),x,repmat(x(:,n),p,ceil(K/2))];
y = zeros(p,n);
if n>1e4 && verbose
    fprintf(1,'Applying median filter...')
    fprintf(1,'%6.2f %%',0)
end
for i=1:n
    ind = i:i+K;
    y(:,i) = median(x(:,ind),2);
    if mod(i,1e4)<1 && verbose
        fprintf(1,repmat('\b',1,8))
        fprintf(1,'%6.2f %%',100*i/n)
    end
end
if n>1e4 && verbose
    fprintf(1,repmat('\b',1,8))
    fprintf(1,[' OK.'])
    fprintf(1,'\n')
end