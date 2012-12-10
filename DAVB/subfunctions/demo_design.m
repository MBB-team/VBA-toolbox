% demo for design optimization

n =1e2;
N = 64;
gridt = 0:1e-2:1;

mo = zeros(length(gridt),1);
so = zeros(length(gridt),1);
c = [eye(2);ones(1,2)];

kernel = exp(-0.5*(-10:.1:10).^2*4e2);
kenerl = kernel./sum(kernel);
figure,plot(kernel,'.')

for i=1:length(gridt)
    opt = zeros(N,1);
    for j=1:N
        X = [];
        X1 = randn(n,1)>gridt(i);
        X1 = conv(double(X1),kernel','same');
%         X1 = (X1-mean(X1))./std(X1);
        X2 = randn(n,1)>gridt(i);
        X2 = conv(double(X2),kernel','same');
%         X2 = (X2-mean(X2))./std(X2);
        X = [X1,X2,ones(n,1)];
        Cov = pinv(X'*X);
        opt(j) = trace(c'*Cov*c);
    end
    mo(i) = mean(opt);
    so(i) = std(opt)./sqrt(N);
end


figure,errorbar(mo,so)


