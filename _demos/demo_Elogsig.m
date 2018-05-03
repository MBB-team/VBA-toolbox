% demo for Elogsig

gridM = -10:10;
gridV = 0:10:100;
Nmc = 1e4;
for i=1:length(gridM)
    m = gridM(i);
    for j=1:length(gridV)
        V = gridV(j);
        X = m + sqrt(V).*randn(Nmc,1);
        lsx1 = log(sig(X));
        lsx2 = log(1-sig(X));
        Els1(i,j) = mean(lsx1);
        Els2(i,j) = mean(lsx2);
        Els01(i,j) = Elogsig(m,V);
        Els02(i,j) = Elogsig(-m,V);
    end
end

figure,plot(vec(exp(Els1)),vec(exp(Els01)),'.')
figure,plot(vec(exp(Els2)),vec(exp(Els02)),'.')