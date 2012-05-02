function [ha] = imageMat(A)

% display VB score for the missing region problem

hf = figure;
ha = axes('parent',hf);
imagesc(A,'parent',ha);
colorbar('peer',ha)

xlabel(ha,'sampling rate')
ylabel(ha,'SNR')

xl = {'0.5 sec','1 sec','2 sec'};
yl = {'1','10','100'};

set(ha,...
    'xtick',[1 2 3],'xticklabel',xl,...
    'ytick',[1 2 3],'yticklabel',yl)