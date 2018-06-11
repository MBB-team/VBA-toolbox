function [p,out,hf] = doROC(xp,xn,displayWin)
% implements a ROC analysis on samples from alternative hypotheses
% function [p,out] = doROC(xp,xn)
% IN:
%   - xp: values under the 'positive' hypothesis (H1)
%   - xn: values under the 'negative' hypothesis (H0)
% OUT:
%   - p: area under the ROC curve (probability of discriminatng between the
%   two hypotheses)
%   - out: structure with the following fields:
%       .TP: true positive rate for each threshold
%       .TN: true negative rate ...
%       .FP: false positive rate ...
%       .FN: false negative rate ...
%       .threshold: 5% FPR threshold
%       .power: TPR at 5% FPR threshold
%       .t_star: disambiguation threshold. This is the value for which the
%       total error rate (P(e)=FPR+FNR) is minimal.
%   - hf: the figure handle

try;displayWin;catch;displayWin=1;end

ux = unique([xp(:);xn(:)]);
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

% area under the ROC curve
p = -trapz(1-out.TN,out.TP);

[hp,gp] = VBA_empiricalDensity(xp(:),1);
[hn,gn] = VBA_empiricalDensity(xn(:),1);

% find FPR=0.05 threshold and power
d = (0.05-out.FP).^2;
[mfp,mi] = min(d);
out.threshold = ux(mi);
out.power = out.TP(mi);

% find disambiguation threshold
r = out.FP+out.FN;
[mfp,mi] = min(r);
out.t_star = ux(mi);

if ~displayWin
    return
end

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
title(ha(2),['ROC curve: p = ',num2str(p),' , [FPR=0.05]: power=',num2str(out.power),' ,threshold=',num2str(out.threshold)])

ha(3) = subplot(3,1,3,'parent',hf,'xlim',[min(ux),max(ux)],'ylim',[min(r),max(r)],'nextplot','add');
plot(ha(3),ux,r,'b.')
plot(ha(3),[out.t_star,out.t_star],[min(r),max(r)],'k--')
xlabel(ha(3),'threshold values')
ylabel(ha(3),'P(e) = FPR+FNR')
title(ha(3),['disambiguation threshold=',num2str(out.t_star)])


