function h = plotElipse(center,axesLength,hParent)

% This function plots 2D-ellipses in given axes

try
    haf = hParent;
    isax = isequal(get(haf,'type'),'axes');
    if ~isax
        hf = figure('color',ones(1,3));
        haf = axes('parent',hf,'nextplot','add');
        box(haf,'on')
    end
catch
    hf = figure('color',ones(1,3));
    haf = axes('parent',hf,'nextplot','add');
    box(haf,'on')
end

h = rectangle(...
    'Position',[center(:)'-0.5*axesLength(:)',axesLength(:)'],...
    'Curvature',[1,1],...
    'facecolor',0.8*ones(1,3));
