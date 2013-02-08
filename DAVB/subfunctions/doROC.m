function [p,out] = doROC(xp,xn)

ux = unique([xp(:),xn(:)]);
Np = length(xp);
Nn = length(xn);
n = length(ux);
out.FP = zeros(n,1);
out.TP = zeros(n,1);
out.FN = zeros(n,1);
out.TN = zeros(n,1);

for i=1:n
    t = ux(i);
    out.TP(i) = length(find(xp>=t))./Np;
    out.FN(i) = length(find(xp<t))./Np;
    out.FP(i) = length(find(xn>=t))./Nn;
    out.TN(i) = length(find(xn<t))./Nn;
end

[hp,gp] = empiricalHist(xp(:),1);
[hn,gn] = empiricalHist(xn(:),1);



hf = figure('color',[1 1 1]);

ha = subplot(2,1,1,'parent',hf,'nextplot','add');
plot(ha,gp,hp,'g')
plot(ha,gn,hn,'r')

ha = subplot(2,1,2,'parent',hf,'xlim',[0,1],'ylim',[0,1],'nextplot','add');
plot(ha,[0,1],[0,1],'r--')
plot(ha,1-out.TN,out.TP,'k.')
xlabel(ha,'1-TN rate')
ylabel(ha,'TP rate')
p = -trapz(1-out.TN,out.TP);
title(ha,['ROC curve: p = ',num2str(p)])



