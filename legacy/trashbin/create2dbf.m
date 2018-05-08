function X = create2dbf(g1,g2,N)

n1 = length(g1);
n2 = length(g2);

X1 = zeros(n1,N,2);
for i = 1:N
    X1(:,i,1) = cos((i-1).*pi.*[0:(n1-1)]./(n1-1));
    X1(:,i,2) = sin((i-1).*pi.*[0:(n1-1)]./(n1-1));
end
X2 = zeros(n2,N,2);
for i = 1:N
    X2(:,i,1) = cos((i-1).*pi.*[0:(n2-1)]./(n2-1));
    X2(:,i,2) = sin((i-1).*pi.*[0:(n2-1)]./(n2-1));
end

X = zeros(n1,n2,N^2);
ij = 1;
for i=1:N
    for j=1:N
        for i1=1:2
            for i2=1:2
                X(:,:,ij) = repmat(X1(:,i,i1),1,n2).*repmat(X2(:,j,i2)',n1,1);
                ij = ij +1;
            end
        end
    end
end
