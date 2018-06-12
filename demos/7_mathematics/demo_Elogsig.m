% demo for Elogsig

gridM = -10:10;
gridV = 0:10:100;
Nmc = 1e4;
for i=1:length(gridM)
    m = gridM(i);
    for j=1:length(gridV)
        V = gridV(j);
        X = m + sqrt(V).*randn(Nmc,1);
        lsx1 = log(VBA_sigmoid(X));
        lsx2 = log(1-VBA_sigmoid(X));
        Els1(i,j) = mean(lsx1);
        Els2(i,j) = mean(lsx2);
        Els01(i,j) = VBA_Elogsig(m,V);
        Els02(i,j) = VBA_Elogsig(-m,V);
    end
end

figure,plot(VBA_vec(exp(Els1)),VBA_vec(exp(Els01)),'.')
figure,plot(VBA_vec(exp(Els2)),VBA_vec(exp(Els02)),'.')