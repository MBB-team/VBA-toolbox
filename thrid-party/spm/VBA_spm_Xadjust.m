function [X0,isYout] = VBA_spm_Xadjust(SPMfile,VOIfile,varthresh)
% gets the effects of no interest, wrt to which the data has been adjusted
% [X0,isYout] = spm_Xadjust(SPMfile,VOIfile,varthresh)
% IN:
%   - SPMfile: name of the SPM file
%   - VOIfile: name of the VOI file
%   - varthresh: threshold for PCA on X0 (fraction of explained variance).
%   Default is 0.95.
% OUT:
%   - X0: confounds matrix
%   - isYout: vector of indices of scans that were effectively removed from
%   the GLM using scan-nulling regressors included in the original
%   confounds matrix.

try, varthresh; catch, varthresh = 0.95; end

load(SPMfile)
load(VOIfile)

if isequal(xY.Ic,0)
    X0 = [];
    isYout = [];
    return
end
Fc = SPM.xCon(xY.Ic);
ind = sum(abs(Fc.c),2)==0;
X0 = SPM.xX.X(:,ind);

% remove scan-nulling regressors
n0 = size(X0,2);
remove = [];
isYout = [];
for i=1:n0
    X0i = X0(:,i) - X0(1,i);
    i1 = find(X0i~=0);
    tmp = zeros(size(X0,1),1);
    tmp(i1) = 1;
    if isequal(X0i,tmp)
        remove = [remove;i];
        isYout = [isYout,i1];
    end
end
X0 = X0(:,setdiff(1:size(X0,2),remove));
X0 = bsxfun(@minus,X0,mean(X0,1));
X0 = bsxfun(@rdivide,X0,std(X0,[],1));
if varthresh < 1
    [u,s,v] = svd(X0);
    s2 = diag(s).^2;
    ev = cumsum(s2./sum(s2));
    it = find(ev>=varthresh, 1 );
    X0 = [u(:,1:it),ones(size(X0,1),1)./sqrt(size(X0,1))];
end

return

% the following reproduces the code in spm_regions.m:

% sX = SPM.xX.xKXs;
% X0 = sX.X*(eye(spm_sp('size',sX,2)) - spm_sp('xpx-',sX)*sf_H(Fc,sX));
%
% function H = sf_H(Fc,sX)
% if sf_ver(Fc) > 1,
%    hsqr = sf_Hsqr(Fc,sX);
%    H    = hsqr' * hsqr;
% else
%    H = Fc.c * pinv(Fc.X1o' * Fc.X1o) * Fc.c';
% end
%
% function hsqr = sf_Hsqr(Fc,sX)
% if sf_ver(Fc) > 1,
%    hsqr = spm_sp('ox',spm_sp('set',Fc.X1o.ukX1o))' * spm_sp('cukx',sX);
% else
%    hsqr = spm_sp('ox',spm_sp('set',Fc.X1o))'*spm_sp('x',sX);
% end
%
% function v = sf_ver(Fc)
% if isstruct(Fc.X0), v = 2; else v = 1; end

