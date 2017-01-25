function [all] = VBA_classification_display(all)
% displays the results of VBA_classification
% function [all] = VBA_classification_display(all)
% Note: inputs to VBA_classification_display.m are outputs of
% VBA_classification.m.
% IN:
%   - pv: cross-validation p-value
%   - df: effective degees of freedom of the t- or F- test
%   - all: structure array with fields:
%       .P: pxk matrix of classification weights
%       .posterior: VBA's posterior structure for full-data inversion
%       .out: VBA's out structure for full data-inversion
% OUT:

try
    acc = all.stat.success;
    nok = sum(acc);
    n = length(acc);
    r = all.r;
    p = size(all.in.X,2);
    k = all.in.k;
    pv = all.stat.pv;
    pa = all.stat.pa;
    bpa = all.stat.bpa;
    pBayes = all.stat.pBayes;
catch
    disp('Error: ''all'' structure does not contain appropriate entries!')
    return
end

figname = 'VBA classification';
if all.in.sparse
    figname = [figname,' (sparse inversion)'];
end
figname = [figname,': cross-validation results'];
all.hf = figure('color',[1 1 1],'name',figname,'menubar','none');

% display # observed accurate predictions (and its distribution under H0)
ha = subplot(2,2,1,'parent',all.hf,'nextplot','add');
[pH0] = VBA_binomial(0:n,n,r);
bar(0:nok-1,pH0(1:nok),'parent',ha,'facecolor',0.8*[1 1 1]);
bar(nok:n,pH0(nok+1:n+1),'parent',ha,'facecolor',0.8*[1 0.75 0.75]);
plot(nok*[1 1],get(ha,'ylim'),'r')
set(ha,'xlim',[-0.5,n+1.5])
xlabel(ha,'x = number of accurate predictions')
ylabel(ha,'p(x|H_0)')
title(ha,'classical inference')

% display classification weights
ha = subplot(2,2,2,'parent',all.hf,'nextplot','add');
plot(ha,all.P,'marker','.')
set(ha,'xlim',[1,size(all.P,1)])
xlabel(ha,'weights'' dimensions')
title(ha,'weights'' estimates')

% % display correlation among k-folds
% C = corrcoef(all.P);
% ha = subplot(2,2,3,'parent',all.hf,'nextplot','add');
% imagesc(C,'parent',ha);
% colorbar('peer',ha)
% set(ha,'clim',[-1,1])
% axis(ha,'equal')
% axis(ha,'tight')
% title(ha,'estimation stability (among k-folds)')
% xlabel(ha,'k-folds')
% ylabel(ha,'k-folds')

% display bayesian posterior pdf on classification accuracy
ha = subplot(2,2,3,'parent',all.hf,'nextplot','add');
VBA_PPM(nok,n-nok,0.5,'beta',1,ha);
set(ha,'ytick',[])
xlabel(ha,'r = classifier accuracy')
ylabel(ha,'p(r|x)')
title(ha,'Bayesian inference')

% display summary of results
if floor(all.dt./60) == 0
    timeString = [num2str(floor(all.dt)),' sec'];
else
    timeString = [num2str(floor(all.dt./60)),' min'];
end
str{1} = sprintf(['Date: ',datestr(all.date)]);
str{2} = sprintf(['Cross-validation complete (took ~',timeString,')']);
str{3} = sprintf(['Dimensions of the analysis:','\n ',...
    '    - data: n=',num2str(n),'\n ',...
    '    - classification weights: p=',num2str(p),'\n ',...
    '    - # folds: k=',num2str(k),' (training sets overlap = ',num2str((k-2)/(k-1),2),')']);
str{4} = sprintf(['Summary statistics:','\n ',...
    '    - Classical inference: p-val. = ',num2str(pv,3),'\n ',...
    '    - Bayesian inference: exc. prob. = ',num2str(pBayes,3),'\n ',...
    '    - class. acc. = ',num2str(nok),'/',num2str(n),' = ',num2str(pa,3),'\n ',...
    '    - balanced class. acc. = ',num2str(bpa,3)]);
VBA_disp(str(2:4),struct('verbose',1))
disp(' ')
uicontrol('parent',all.hf,'style','text','units','normalized','position',[0.52,0.05,0.5,0.4],'backgroundcolor',[1,1,1],'HorizontalAlignment','left','fontsize',9,'string',str);

try,getSubplots;end


