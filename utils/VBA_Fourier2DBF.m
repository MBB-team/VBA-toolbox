function [X,X1,X2] = VBA_Fourier2DBF(L1,L2,N,dct)
% creates 2D Fourier basis functions
% function X = create2dbf(g1,g2,N)
% This function evaluates Fourier basis function sets, evaluated on a 2D
% domain defined by 2 grids g1 and g2. Let z(x) be a 1D-field of length L.
% It can be decomposition using the Following Fourier series:
% z(x)=a0+sum_n[an*cos(2*n*pi*x/L)+bn*sin(2*n*pi*x/L)]
% where an and bn are the coefficient of the series.
% This decomposition can be generalized to a 2D-field z(x,y) of size
% [L1,L2], as follows:
% z(x,y)=a0+sum_n[sum_m[anm*fn(x,y)+bnm*gnm(x,y)+cnm*hnm(x,y)+inm*f4(x,y)]]
% where the Basis functions are defined as follows:
% fnm(x,y) = cos(2*n*pi*x/L1)*cos(2*m*pi*y/L2)
% gnm(x,y) = cos(2*n*pi*x/L1)*sin(2*m*pi*y/L2)
% hnm(x,y) = sin(2*n*pi*x/L1)*cos(2*m*pi*y/L2)
% inm(x,y) = sin(2*n*pi*x/L1)*sin(2*m*pi*y/L2)
% An almost equivalent alternative is DCT (Discrete Cosinus Transform),
% which only uses cos basis functions, but with finer (double) frequency
% sampling in the set.
% IN:
%   - L1: grid for 1st dimension of the domain
%   - L2: grid for 2nd dimension of the domain
%   - N: order of the Fourier series
%   - dct: flag for DCT (dct=1) or Fourier (dct=0, default) basis set.
% OUT:
%   - X: L1xL2x(4*N*(N+1)+1) array, where X(:,:,i) is the i^th basis
%   function.
%   - X1: L1x(2*N)x2 array, where X1(:,i,1) is cos(2*n*pi*x/L1) and
%   X1(:,i,2) is sin(2*n*pi*x/L1).
%   - X2: [id for 2nd dim]

try,dct;catch,dct=0;end


% 1D cos- and sin- basis function sets
X1 = zeros(L1,N,2);
X2 = zeros(L2,N,2);
X = zeros(L1,L2,4*N*(N+1)+1);
X(:,:,1) = 1; % constant term
ij = 2;
for i = 1:N
    if dct
        X1(:,i,1) = cos((2*i-1).*pi.*[0:(L1-1)]./(L1-1));
    else
        X1(:,i,1) = sin(2*i.*pi.*[0:(L1-1)]./(L1-1));
    end
    X1(:,i,2) = cos(2*i.*pi.*[0:(L1-1)]./(L1-1));
    if dct
        X2(:,i,1) = cos((2*i-1).*pi.*[0:(L2-1)]./(L2-1));
    else
        X2(:,i,1) = sin(2*i.*pi.*[0:(L2-1)]./(L2-1));
    end
    X2(:,i,2) = cos(2*i.*pi.*[0:(L2-1)]./(L2-1));
    X(:,:,ij) = repmat(X1(:,i,1),1,L2);
    X(:,:,ij+1) = repmat(X1(:,i,2),1,L2);
    X(:,:,ij+2) = repmat(X2(:,i,1)',L1,1);
    X(:,:,ij+3) = repmat(X2(:,i,2)',L1,1);
    ij = ij +4;
end

% 2D Fourier basis function set
for i=1:N
    for j=1:N
        for i1=1:2
            for i2=1:2
                X(:,:,ij) = repmat(X1(:,i,i1),1,L2).*repmat(X2(:,j,i2)',L1,1);
                ij = ij +1;
            end
        end
    end
end

