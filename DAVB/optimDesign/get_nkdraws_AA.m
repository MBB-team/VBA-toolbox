function [C,dt] = get_nkdraws_AA(k,n,verbose)
% gets all n-draws (with replacement) from a k-urn.
% [C] = get_nkdraws(k,n,verbose)
% This function construct all the multivariate draws of size n that can be
% obtained (with replacement) from an urn containing k marbles.
% A couple of examples is worth a long description:
%
% get_nkdraws(2,3)
%
% ans =
%
%      1     1     1     1     2     2     2     2
%      1     1     2     2     1     1     2     2
%      1     2     1     2     1     2     1     2
%
% get_nkdraws(3,2)
%
% ans =
%
%      1     1     1     2     2     2     3     3     3
%      1     2     3     1     2     3     1     2     3
%
% IN:
%   - k: the cardinality of the set, whose permutations are replicated over
%   n dimensions
%   - n: the dimension of the vector, whose entries are elements of the set
%   {1,2,...,k}.
% OUT:
%   - C: a nXn^p matrix, whose columns are the n-draws

try; verbose; catch verbose = 0; end
tic

np = k^n;
if verbose
    fprintf(1,['Getting the ',num2str(np),' ',...
        num2str(n),'-draws from ',num2str(k),'-urn...'])
end
C = zeros(n,np); % pre-allocate k-permutations
for ind = 1:n
    for j = 1:k^(ind-1)
        nbch = np/(k^(ind-1));
        C(ind,(1+(j-1)*nbch):j*nbch) = ceil((1:nbch)/(nbch/k));
    end
end

dt = toc;
if verbose
    fprintf(1,[' OK. (took ',num2str(dt),' seconds)'] )
    fprintf(1,'\n')
end

return
