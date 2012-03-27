function [h] = displayOptimDesign(muy,Vy,u,dim,DJS,b)
% plots the prior predictive densities under different designs
%
% function [h] = displayOptimDesign(muy,Vy,u,dim,DJS,b)
% IN:
%   - muy,Vy: ndxnm cell array of the 1st- and 2nd-order moments of the
%   prior predictive density under the nm models and nd designs
%   - u: ndx1 cell array of different designed DCM inputs
%   - dim: ndxnm cell array of dim structures
%   - DJS: ndx1 vector of Jensen-Shannon divergences for each design
%   - b: ndx1 vector of ensuing lower bound on the probability of
%   classification error
% OUT:
%   - h: handles structure
%------------------------------------------------------------
% Copyright (C) 2012 Jean Daunizeau / License GNU GPL v2
%------------------------------------------------------------

pos0 = get(0,'screenSize');
pos = [0.45*pos0(3),0.05*pos0(4),0.5*pos0(3),0.85*pos0(4)];
h.f = figure(...
    'position',pos,...
    'color',[1 1 1],...
    'name','Design optimization tool');

[nd,nm] = size(muy);
ny = dim{1}.p;
nu = size(u{1},1);

% first plot designs
for i= 1:nd
    h.au(i) = subplot(nm+1,nd,i,...
        'parent',h.f);
    imagesc(u{i},'parent',h.au(i))
    title(h.au(i),...
        {['design ',num2str(i)],...
        [' (DJS=',num2str(DJS(i)),...
        ' ; p=',num2str(b(i)/2),')']})
    xlabel(h.au(i),'time samples')
    ylabel(h.au(i),'inputs')
    set(h.au(i),'ytick',1:nu);
end


for j=1:nm %loops over models

    for i=1:nd % loops over experimental designs
        
        h.av(i,j) = subplot(nm+1,nd,j*nd+i,'parent',h.f);
        imagesc(reshape(diag(Vy{i,j}),ny,[]),'parent',h.av(i,j))
        title(h.av(i,j),...
            ['design ',num2str(i),...
            ' - model ',num2str(j)])
        xlabel(h.av(i,j),'time samples')
        ylabel(h.av(i,j),'channel variance')
        set(h.av(i,j),'ytick',1:ny);
        
    end
    
end

try, getSubplots; end
    
    
    
    
