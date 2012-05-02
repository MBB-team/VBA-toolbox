function [C,dt] = get_nkdraws(k,n,verbose)
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
%   - C: a nXk^n matrix, whose columns are the n-draws

try; verbose; catch verbose = 0; end
tic

strFOR      = [];
strEND      = ['if C(n,it)~=0;it=it+1;end;'];
strIND      = [];
for i = 1:n             % loop on the set elements
    strFOR  = [strFOR,'for i',num2str(i),'=1:k;'];
    strIND  = [strIND,';i',num2str(i)];
    strEND  = [strEND,'end;'];
end

np = k^n;
if verbose
    fprintf(1,['Getting the ',num2str(np),' ',...
        num2str(n),'-draws from ',num2str(k),'-urn...'])
end
C = zeros(n,np); % pre-allocate k-permutations
it = 1; % initialize column-index within the loop
eval([strFOR,'C(:,it)=[',strIND,'];',strEND]) % loop

dt = toc;
if verbose
    fprintf(1,[' OK. (took ',num2str(dt),' seconds)'] )
    fprintf(1,'\n')
end

return


