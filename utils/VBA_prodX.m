function [gx, dgdx, d2gdx2] = VBA_prodX (X1, X2)
% // VBA toolbox //////////////////////////////////////////////////////////
%
% [gx, dgdx, d2gdx2] = VBA_prodX (X1, X2)
% calculates the 1st and 2nd order derivatives of the dot product X1' * X2
% of two vectors X1 and X2
%
% IN:
%   - X1,X2: two vectors of same size (one of them can be a scalar)
% OUT:
%   - gx: the product of the two vectors (X1'*X2)
%   - dgdx: the 1st order derivative of the product
%   - d2gdx2: the 2nd order derivative of the product
%
% /////////////////////////////////////////////////////////////////////////

if ~isequal(size(X1),size(X2))
    if ~isequal(numel(X1),1) && ~isequal(numel(X2),1)
        disp('Error: entries must have same size or one of them must be scalar.')
        gx = [];
        dgdx = [];
        d2gdx2 = [];
        return
    end
end

gx = sum(X1'*X2);

n1 = size(X1,1);
n2 = size(X2,1);

if isequal(n1,n2)
    dgdx = [X2;X1];
    d2gdx2 = [  zeros(n2,n1)    eye(n2)
                eye(n1)         zeros(n1,n2)  ];
elseif isequal(n1,1)
    dgdx = [sum(X2);X1.*ones(n2,1)];
    d2gdx2 = [  0           ones(1,n2)
                ones(n2,1)  zeros(n2,n2)  ];
else
    dgdx = [X2.*ones(n1,1);sum(X1)];
    d2gdx2 = [  zeros(n1,n1)    ones(n1,1)
                ones(1,n1)  	0           ];
end




