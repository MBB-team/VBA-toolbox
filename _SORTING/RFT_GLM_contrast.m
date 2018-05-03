function [stat,out] = RFT_GLM_contrast(X,y,c,type,u,verbose)
% applies RFT to GLM-based contrast inference
% function [stat,out] = RFT_GLM_contrast(X,y,c,type,u,verbose)
%
% In brief, this function uses RFT to test linear mixtures of effects,
% under the general linear model (GLM):
%   y = X*beta + e
% where y is the data, e is some (Gaussian) random noise, X is the design
% matrix and beta are unknown GLM parameters.
% Tests are specified in terms of a contrast matrix c, such that the
% p-value computes the probability of the test score under the null H0.
% More precisely:
%  - 't' tests are one-sided: H0={c'*beta>0}
%  - 'F' tests are two-sided: H0={c(:,1)'beta=0 AND ... c(:,m)'*beta=0}
% For example, testing for the ith parameter corresponds to an all-zero
% contrast vector, except on its ith entry.
% Note that here, y is repeatedly sampled over a lattice, i.e. it is a
% matrix of size nxL, where L is the number of positions on the lattice.
% RFT is then applied to derive p-values on topological features of the
% statistical random field sampled at each position on the lattice.
%
% IN:
%   - X: nXk design matrix
%   - y: nXL data matrix
%   - c: kXm contrast matrix (default is eye(k) -> omnibus F-test)
%   - type: flag for t- or F- test. Can be set to 't' (default) or 'F'
%   - u: set-inducing threshold (for cluster inference)
%   - verbose: flag for displaying results (default is 0)
% OUT:
%   - stat: vector containing the t- or F-values
% (depending on the test type selected in the input) corresponding to each point in Y
%   - out: RFT output structure (see RFT_main.m).
% 
% See also RFT_main.m

% fill in default I/O
out = [];
[n,L] = size(y);
[n0,k] = size(X);
[k0,m] = size(c);
try;c;catch;c=eye(k);type='F';end
try;type;catch;type='t';end
try;u;catch;u=[];end
try;verbose;catch;verbose=0;end

% check basic numerical requirements
try
    if VBA_isWeird (y)
        disp('Error: data contains weird values!')
        return
    end
end
if ~isequal(n,n0)
    disp('Error: design matrix has to have as many rows as the data matrix!')
    return
end
if ~isequal(k,k0)
    disp('Error: contrasts have to have as many rows as there are columns in the design matrix!')
    return
end


% Estimate GLM parameters
if verbose
    fprintf(1,'Computing parameter covariance matrix...')
end
C = X'*X;
iC = pinv(C);
if verbose
    fprintf(1,' OK.')
    fprintf(1,'\n')
end
iX = iC*X'; % pseudo-inverse of design matrix X
b = iX*y;
if verbose
    fprintf(1,'Computing projection matrices...')
end
P = X*iX;
R = speye(n) - P;
if verbose
    fprintf(1,' OK.')
    fprintf(1,'\n')
end
yhat = X*b;
trR = n - trace(P);
df = trR.^2./sum(sum(R.^2,1));

% estiamte statistical field
switch type
    case 't'
        if m~=1
            disp('Error: cannot have contrast matrices for ''t''-tests!')
            return
        end
        vhat = sum((y-yhat).^2,1)./trR;
        V = vhat.*(c'*iC*c);
        stat = (c'*b)./sqrt(V);
        stat = stat';
    case 'F'
        if verbose
            fprintf(1,'Computing contrast null-space projectors...')
        end
        ic = pinv(c'*c)*c';
        c0 = speye(size(c,1)) - c*ic;
        X0 = X*c0;
        R0 = speye(n) - X0*pinv(X0'*X0)*X0';
        M = R0 - R;
        trM = trace(M);
        df = [trM.^2./sum(sum(M.^2,1)),df];
        if verbose
            fprintf(1,' OK.')
            fprintf(1,'\n')
        end
        if verbose
            fprintf(1,'Computing F-statistics...')
            if L>1
                fprintf(1,'%6.2f %%',0)
            end
        end
        stat = zeros(L,1);
        for i=1:L
            stat(i) = ((yhat(:,i)'*M*yhat(:,i))./(y(:,i)'*R*y(:,i))).*(trR./trM);
            if verbose && L>1
                fprintf(1,repmat('\b',1,8))
                fprintf(1,'%6.2f %%',100*i/L)
            end
        end
        if verbose
            if L>1
                fprintf(1,repmat('\b',1,8))
            end
            fprintf(1,[' OK.'])
            fprintf(1,'\n')
        end
    otherwise
        disp('Error: this function only supports t- and F- tests!')
        return
end

% call RFT
% Z = RFT_Gtf(stat,df);
options.R = (y-yhat)';
options.type = type;
options.dof = df;
if ~isempty(u)
    options.u = u;
end
[out] = RFT_main(stat,options,verbose);



