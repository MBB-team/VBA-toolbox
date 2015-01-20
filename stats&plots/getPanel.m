function hfig = getPanel(hfig)
    try
        get(hfig,'BackgroundColor');
    catch % in 2014b gco parent is the figure not the panel...
        c = get(hfig,'Children');
        isPanel = cell2mat(arrayfun(@(x) strcmp(class(x),'matlab.ui.container.Panel'),c,'uniformOutput',false));
        hfig = c(isPanel);
    end
