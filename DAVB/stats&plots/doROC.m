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

% find FPR=0.05 threshold and power
d = (0.05-out.FP).^2;
[mfp,mi] = min(d);
out.threshold = ux(mi);
out.power = out.TP(mi);

% find reversal threshold
r = out.FP+out.FN;
[mfp,mi] = min(r);
out.t_star = ux(mi);

hf = figure('color',[1 1 1]);

ha = subplot(3,1,1,'parent',hf,'nextplot','add');
plot(ha,gp,hp,'g')
plot(ha,gn,hn,'r')
plot(ha,out.t_star*[1,1],[0,max([hp(:);hn(:)])],'k--')
title(ha,'empirical histograms')
legend({'H1','H0'})

ha(2) = subplot(3,1,2,'parent',hf,'xlim',[0,1],'ylim',[0,1],'nextplot','add');
plot(ha(2),[0,1],[0,1],'r--')
plot(ha(2),1-out.TN,out.TP,'b.')
plot(ha(2),[0.05,0.05],[0,out.power],'r')
plot(ha(2),[0,0.05],[out.power,out.power],'r')
plot(ha(2),0.05,out.power,'r.')
xlabel(ha(2),'FPR (1-specificity)')
ylabel(ha(2),'TPR (sensitivity)')
p = -trapz(1-out.TN,out.TP);
title(ha(2),['ROC curve: p = ',num2str(p),' , [FPR=0.05]: power=',num2str(out.power),' ,threshold=',num2str(out.threshold)])

ha(3) = subplot(3,1,3,'parent',hf,'xlim',[min(ux),max(ux)],'ylim',[min(r),max(r)],'nextplot','add');
plot(ha(3),ux,r,'b.')
plot(ha(3),[out.t_star,out.t_star],[min(r),max(r)],'k--')
xlabel(ha(3),'threshold values')
ylabel(ha(3),'P(e) = FPR+FNR')
title(ha(3),['reversal threshold=',num2str(out.t_star)])


