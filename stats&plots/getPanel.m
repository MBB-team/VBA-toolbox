function hpanel = getPanel(hfig)

    c = get(hfig,'Children');
    isPanel = cell2mat(arrayfun(@(x) strcmp(get(x,'Type'),'uipanel'),c,'uniformOutput',false));
    hpanel = c(isPanel);

