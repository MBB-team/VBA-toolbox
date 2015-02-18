function hfig = getPanel(hfig)
    if verLessThan('matlab','8.4.0')
        return
    else % in 2014b gco parent is the figure not the panel...
%         hpanel = findobj(hfig,'Type','Panel');
        c = get(hfig,'Children');
        isPanel = cell2mat(arrayfun(@(x) strcmp(class(x),'matlab.ui.container.Panel'),c,'uniformOutput',false));
        hpanel = c(isPanel);
        if ~isempty(hpanel)
            hfig = hpanel ;
        end
        
    end
