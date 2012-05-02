% OTO: script for OTO results

% first get real association over time
load('D:\MatlabWork\routinetheque\Learning\OTO\u.mat')

% then get perceived association
x = posterior.muX(2,:);

p11 = sigm(p1(indC1)./100,struct('INV',1));
% p21 = p1(indC2);

vals = unique(p11);
n = length(vals);

yy = [];
xx = [];

for i=1:n
    
    indi = find(p11==vals(i));
    tmp = x(indi);
    yy = [yy;tmp(:)];
    xx = [xx,p11(indi)];
    m(i) = mean(tmp);
    s(i) = std(tmp);
    
end

miy = min([xx(:);yy(:)]);
may = max([xx(:);yy(:)]);


figure
plot(xx,yy,'b.')
hold on
plot([miy,may],[miy,may],'r')
errorbar(vals,m,s,'r.')
title('controlled vs perceived association')
xlabel('controlled')
ylabel('perceived')
grid on
axis tight