function colors = getColors(n)
hf = figure('visible','off');
ha = axes('parent',hf);
colors = get(ha,'colororder');
colors = repmat(colors,ceil(n/size(colors,1)),1);
colors = colors(1:n,:);
delete(hf)