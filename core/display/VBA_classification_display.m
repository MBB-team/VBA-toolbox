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
    Eg = all.Eg;
    Vg = all.Vg;
catch
    disp('Error: ''all'' structure does not contain appropriate entries!')
    return
end

pos0 = get(0,'screenSize');
pos = [0.51*pos0(3),0.05*pos0(4),0.45*pos0(3),0.9*pos0(4)];
figname = 'VBA classification';
if all.in.sparse
    figname = [figname,' (sparse inversion)'];
end
figname = [figname,': cross-validation results'];
all.handles.hf = figure('position',pos,'color',[1 1 1],'name',figname,'menubar','none','tag','VBA_classification');

% display classifier prediction
all.handles.ha(1) = subplot(4,2,1,'parent',all.handles.hf,'nextplot','add');
plot(all.handles.ha(1),Eg,'k')
% plotUncertainTimeSeries(Eg',Vg',[],all.handles.ha(1),[],'k');
plot(all.handles.ha(1),[0.5,length(acc)+0.5],[0.5,0.5],'k--')
plot(all.handles.ha(1),find(acc==1),Eg(acc==1),'g.')
plot(all.handles.ha(1),find(acc==0),Eg(acc==0),'r.')
set(all.handles.ha(1),'xlim',[0.5,length(acc)+0.5])
xlabel(all.handles.ha(1),'data samples')
ylabel(all.handles.ha(1),'test prediction E[y]')
VBA_title(all.handles.ha(1),'classifier test prediction')

% display classification accuracy
all.handles.ha(2) = subplot(4,2,2,'parent',all.handles.hf,'nextplot','add');
[stacky,stdy,gridg,stdg] = VBA_Bin2Cont(Eg,all.in.y);
plot(all.handles.ha(2),[0,1],[0,1],'r')
gridp = 0:1e-2:1;
plot(all.handles.ha(2),gridp,gridp+sqrt(gridp.*(1-gridp)),'r--')
plot(all.handles.ha(2),gridp,gridp-sqrt(gridp.*(1-gridp)),'r--')
errorbar(gridg,stacky,stdy,'k.','parent',all.handles.ha(2))
grid(all.handles.ha(2),'on')
axis(all.handles.ha(2),'tight')
xlabel(all.handles.ha(2),'test prediction E[y]')
ylabel(all.handles.ha(2),'data average E[y]')
VBA_title(all.handles.ha(2),'test prediction accuracy')

% display classifier doubt
iok = find(acc==1);
mvok = mean(Vg(iok));
vvok = var(Vg(iok));
inok = find(acc==0);
mvnok = mean(Vg(inok));
vvnok = var(Vg(inok));
all.handles.ha(3) = subplot(4,2,3,'parent',all.handles.hf,'nextplot','add');
[haf,hf,hp] = plotUncertainTimeSeries([mvok;mvnok],[vvok;vvnok],[],all.handles.ha(3),[],0.8*[1 1 1]);
set(hp,'barwidth',0.4)
set(hf,'color','k')
xlabel(all.handles.ha(3),'classification outcomes')
ylabel(all.handles.ha(3),'prediction variance V[y]')
VBA_title(all.handles.ha(3),'test prediction doubt')
set(all.handles.ha(3),'ylim',[0,0.25],'xlim',[0.5,2.5],'xtick',1:2,'xticklabel',{'correct','wrong'})

% display classification weights
all.handles.ha(4) = subplot(4,2,4,'parent',all.handles.hf,'nextplot','add');
plot(all.handles.ha(4),[1,size(all.P,1)],[0,0],'k--')
plot(all.handles.ha(4),all.P,'marker','.')
set(all.handles.ha(4),'xlim',[0.5,size(all.P,1)+0.5],'xtick',1:size(all.P,1))
xlabel(all.handles.ha(4),'weights'' dimensions')
VBA_title(all.handles.ha(4),'weights'' estimates')

% display # observed accurate predictions (and its distribution under H0)
all.handles.ha(5) = subplot(4,2,5,'parent',all.handles.hf,'nextplot','add');
[pH0] = VBA_binomial(0:n,n,r);
bar(0:nok-1,pH0(1:nok),'parent',all.handles.ha(5),'facecolor',0.8*[1 1 1]);
bar(nok:n,pH0(nok+1:n+1),'parent',all.handles.ha(5),'facecolor',0.8*[1 0.75 0.75]);
plot(all.handles.ha(5),nok*[1 1],get(all.handles.ha(5),'ylim'),'r')
set(all.handles.ha(5),'xlim',[-0.5,n+1.5])
xlabel(all.handles.ha(5),'x = number of accurate predictions')
ylabel(all.handles.ha(5),'p(x|H_0)')
VBA_title(all.handles.ha(5),'classical inference')

% display bayesian posterior pdf on classification accuracy
all.handles.ha(6) = subplot(4,2,6,'parent',all.handles.hf,'nextplot','add');
VBA_PPM(nok+all.out.options.n0,n-nok+all.out.options.n0,r,'beta',1,all.handles.ha(6));
plot(all.handles.ha(6),[r r],get(all.handles.ha(6),'ylim'),'r')
set(all.handles.ha(6),'ytick',[])
xlabel(all.handles.ha(6),'r = classifier accuracy')
ylabel(all.handles.ha(6),'p(r|x)')
VBA_title(all.handles.ha(6),'Bayesian inference')

% add uicontrol for eyeballing VBA inversion results
all.handles.uic = uicontrol('parent',all.handles.hf,'style','pushbutton',...
    'units','normalized','position',[0.4,0.96,0.2,0.02],'backgroundcolor',...
    .8*[1,1,1],'string','eyeball VBA inversion?',...
    'callback',@diagnoseVBA);


% display summary of results
if floor(all.dt./60) == 0
    timeString = [num2str(floor(all.dt)),' sec'];
else
    timeString = [num2str(floor(all.dt./60)),' min'];
end
str{1} = sprintf(['Date: ',datestr(all.date)]);
str{2} = sprintf(['Cross-validation complete (took ~',timeString,')']);
str{3} = sprintf(['Dimensions of the analysis:','\n ',...
    '    - data: n = ',num2str(n),'\n ',...
    '    - classification weights: p = ',num2str(p),'\n ',...
    '    - folds: k = ',num2str(k),' (training sets overlap = ',num2str((k-2)/(k-1),2),')']);
str{4} = sprintf(['Data set properties:','\n ',...
    '    - class imbalance: P(y=1) = ',num2str(r),'\n ',...
    '    - features'' tolerance: tol = [ ',num2str(GLM_tolerance(all.in.X)',2),' ]']);
str{5} = sprintf(['Summary statistics:','\n ',...
    '    - Classical inference: p-val. = ',num2str(pv,3),'\n ',...
    '    - Bayesian inference: exc. prob. = ',num2str(pBayes,3),'\n ',...
    '    - classification accuracy = ',num2str(nok),'/',num2str(n),' = ',num2str(pa,3),'\n ',...
    '    - balanced classification accuracy = ',num2str(bpa,3)]);

for i=1:length(str)
    str{i} = sprintf([str{i},'\n ']);
end
VBA_disp(str(2:5),struct('verbose',1))
disp(' ')
all.handles.huic = uicontrol('parent',all.handles.hf,'style','text','units',...
    'normalized','position',[0.1,0.02,0.8,0.25],'backgroundcolor',[1,1,1],...
    'HorizontalAlignment','left','fontsize',10,'string',str([1,3:5]));

set(all.handles.hf,'userdata',all)


try, VBA_getSubplots (); end



function diagnoseVBA(ho,e)
ud = get(get(ho,'parent'),'userdata');
posterior = ud.posterior;
out = ud.out;
out.options.figName = 'VBA classification: model inversion given whole-data';
VBA_ReDisplay(posterior,out,1);


