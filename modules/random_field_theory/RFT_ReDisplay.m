function [out] = RFT_ReDisplay(X,out)
% display RFT-1D analysis results
% function [out] = RFT_ReDisplay(X,out)
% When called, this code will display the sampled RF field along with
% summary statistics of the RFT analysis. Right-clicking on the field's
% topological features (i.e. clusters and/or local peaks) will display the
% corresponding RFT summary statistics. 
% A note on the colour code: upcrossing clusters will be displayed in green
% if they pass the (user-specified) set-inducing extent threshold, and in
% red otherwise.
% IN:
%   - X: the sampled RF field
%   - out: the output of RFT_main.m
% OUT:
%   - out: the out structure now contains the handle of the results figure
%   (stored in out.hf).

L = length(X);
peaks = out.peaks;
clusters = out.clusters;
options = out.options;
try,nc=length(clusters.k);catch,nc=0;end
OUTSTR = out.OUTSTR;
STR = out.STR;

switch out.options.type
    case 'norm'
        fieldtype = 'gaussian field';
    case 't'
        fieldtype = ['t field with dof=',num2str(out.options.dof)];
    case 'F'
        fieldtype = ['F field with dof=[',num2str(out.options.dof(1)),',',num2str(out.options.dof(2)),']'];
end
figname = ['1D-RFT analysis on ',fieldtype,' (set-inducing thresholds: X>',num2str(options.u,'%3.2f'),', k>',num2str(options.k),')'];

pos0 = get(0,'screenSize');
pos = [0.51*pos0(3),0.05*pos0(4),0.45*pos0(3),0.85*pos0(4)];
out.hf = figure('position',pos,'color',[1 1 1],'name',figname,'tag','RFT1D');

ha = subplot(3,1,1,'parent',out.hf,'nextplot','add','xlim',[1,L]);
plot(ha,X,'k');
plot(ha,[1,L],[options.u,options.u],'--','color',0.8*[1 1 1])
for i=1:nc
    hcmenu = uicontextmenu;
    uimenu(hcmenu, 'Label',['UPCROSSING CLUSTER (X>',num2str(options.u,'%3.2f'),')']);
    uimenu(hcmenu, 'Label',['location: ',num2str(min(clusters.ind{i})),'<t<',num2str(max(clusters.ind{i}))]);
    uimenu(hcmenu, 'Label',['spatial extent: k=',num2str(clusters.k(i))]);
    uimenu(hcmenu, 'Label',['p=',num2str(clusters.prft(i),'%3.3f'),' [RFT]']);
    if out.clusters.prft(i)<out.options.FWER % clusters.k(i)>=options.k
        plot(ha,clusters.ind{i},X(clusters.ind{i}),'g','uicontextmenu',hcmenu);
    else
        plot(ha,clusters.ind{i},X(clusters.ind{i}),'r','uicontextmenu',hcmenu);
    end
end
for i=1:length(peaks.ind)
    hcmenu = uicontextmenu;
    uimenu(hcmenu, 'Label','LOCAL PEAK');
    uimenu(hcmenu, 'Label',['location: t=',num2str(peaks.ind(i))]);
    uimenu(hcmenu, 'Label',['RF value: X=',num2str(peaks.val(i),'%3.3f')]);
    uimenu(hcmenu, 'Label',['p=',num2str(peaks.prft(i),'%3.3f'),' [RFT]']);
    uimenu(hcmenu, 'Label',['p=',num2str(peaks.punc(i),'%3.3f'),' [unc]']);
    if out.peaks.prft(i)<out.options.FWER
        hp = plot(ha,peaks.ind(i),peaks.val(i),'g.','uicontextmenu',hcmenu);
    else
        hp = plot(ha,peaks.ind(i),peaks.val(i),'r.','uicontextmenu',hcmenu);
    end
end
xlabel(ha,'location on the lattice (t)')
ylabel(ha,'random field value (X(t))')

bgc = [1 1 1];
unit = 'normalized';
set(ha,'units',unit);
pos0 = get(ha,'position');
pos = [0.075 pos0(2)-0.25 0.85 0.18];
ht = uicontrol('Style','text','units',unit,'position',pos,'backgroundcolor',bgc,...
    'string',OUTSTR(1:7),'horizontalAlignment','left','FontSize',11,'parent',out.hf);

pos = [0.075 0.04 0.15 pos(2)-0.06];
str = cat(1,{'set: p (c)';'----------------------'},STR.set);
ht = uicontrol('Style','text','units',unit,'position',pos,'backgroundcolor',bgc,...
    'string',str,'horizontalAlignment','left','FontSize',11,'parent',out.hf);

pos(1) = pos(1) + 0.175;
str = cat(1,{'clusters: p (k)';'----------------------'},transpose(STR.cluster));
ht = uicontrol('Style','text','units',unit,'position',pos,'backgroundcolor',bgc,...
    'string',str,'horizontalAlignment','left','FontSize',11,'parent',out.hf);

pos(1) = pos(1) + 0.175;
str = cat(1,{'peaks: p (X)';'----------------------'},transpose(STR.peak));
ht = uicontrol('Style','text','units',unit,'position',pos,'backgroundcolor',bgc,...
    'string',str,'horizontalAlignment','left','FontSize',11,'parent',out.hf);

pos(1) = pos(1) + 0.175;
str = cat(1,{'unc: p';'----------------------'},transpose(STR.unc));
ht = uicontrol('Style','text','units',unit,'position',pos,'backgroundcolor',bgc,...
    'string',str,'horizontalAlignment','left','FontSize',11,'parent',out.hf);

pos(1) = pos(1) + 0.175;
str = cat(1,{'location: t';'----------------------'},transpose(STR.loc));
ht = uicontrol('Style','text','units',unit,'position',pos,'backgroundcolor',bgc,...
    'string',str,'horizontalAlignment','left','FontSize',11,'parent',out.hf);


try, VBA_getSubplots (); end





