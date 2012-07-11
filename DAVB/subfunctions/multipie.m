function h = multipie(X,gridX,gridY,ha)

try
    ha;
    set(ha,'nextplot','add')
catch
    hf = figure('color',[1 1 1]);
    ha = axes('parent',hf,'nextplot','add');
end

[n,npies] = size(X);
% compute min distance between points on the grid
if npies > 1
    d = zeros(n,n);
    for i=1:npies
        for j=i+1:npies
            d(i,j) = sqrt(diff(gridX([i,j])).^2+diff(gridY([i,j])).^2);
        end
    end
    d = d(:);
    scale = 2.1/min(d(d~=0));
else
    scale = 2.1;
end

npies = size(X,2);
colors = getColors(n);
for i=1:npies
    h = pie(X(:,i),'parent',ha);
    ii = 0;
    for j=1:length(h)
        if isequal(get(h(j),'type'),'patch')
            ii = ii+1;
            set(h(j),...
                'xdata',get(h(j),'xdata')/scale+gridX(i),...
                'ydata',get(h(j),'ydata')/scale+gridY(i),...
                'FaceColor',colors(ii,:),...
                'edgecolor','none')
        else
            delete(h(j))
        end
    end
end
axis(ha,'equal')


