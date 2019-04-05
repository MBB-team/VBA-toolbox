function [haf,hf,hp] = plotUncertainTimeSeries(muX,SX,dTime,hParent,ind,color)
% plots uncertain time series
% function [haf,hf,hp] = plotUncertainTimeSeries(muX,SX,dTime,hParent,ind)
% This function plots uncertain time series, ie time series associated with
% their instantaneous variance.
% IN:
%   - muX: a nxT matrix containing the time series
%   - SX: a nxT matrix containing the variances associated with each entry
%   of the muX matrix
%   - dTime: a 1xT vector containing the time values of the decimation of
%   the times series
%   - hParent: the handle of the parent axes
%   - ind: the indices of the dimensions of the time series to be displayed
%   (useful in high dimensional cases)
%   - color: color of the graphical objects
% OUT:
%   - haf: the current axes handle
%   - hf: the patchs handles (for later error bar corrections)
%   - hp: the bar/line plot handles

if nargin < 6
	color=[];
end

% Get dimensions
n = size(muX,1);
if ~exist('dTime','var') || isempty(dTime)
    dTime = 1:size(muX,2);
end
indEnd = length(dTime);

% Get axes parent handle
try
    haf = hParent;
    isax = isequal(get(haf,'type'),'axes');
    if ~isax
        hf = figure('color',ones(1,3));
        haf = axes('parent',hf);
        box(haf,'off');        
    end
catch
    hf = figure('color',ones(1,3));
    haf = axes('parent',hf);
    box(haf,'off')
    if sum(SX(:)) ~= 0
        noButton = 0;
    end
    % Preset axes limits  
end
    
% Get display indices
if ~exist('ind','var') || isempty(ind)
    ind = [1:n];
else
    n = length(ind);
end
sc = 1;
% Display uncertain time series
if indEnd > 1
    
    hp = findobj(haf,'Tag','putLine');
    if isempty (hp) 
        % Plot first moment
        hp = plot(haf,dTime,muX(ind,1:indEnd)','Tag','putLine');
        set(haf,'nextplot','add')
    else
        for i = 1 : n
        set(hp(i), ...
            'XData', dTime, ...
            'YData', muX(ind(i),1:indEnd)'); 
        end
    end
    % Add confidence intervals
    hf = findobj(haf,'Tag','putFill');
    if isempty(hf)
        for i = 1 : n
            yp = [muX(ind(i),1:indEnd)+sc*sqrt(SX(ind(i),1:indEnd)),fliplr(muX(ind(i),1:indEnd)-sc*sqrt(SX(ind(i),1:indEnd)))];
            yp(isnan(yp)) = VBA_nanmean(yp);
            xp = [dTime,fliplr(dTime)];
            col = get(hp(i),'color');           
            hf(i) = fill(xp,yp,'r','parent',haf,'facecolor',col,'edgealpha',0,'facealpha',0.25,'Tag','putFill');
        end
    else
        for i = 1 : n
            yp = [muX(ind(i),1:indEnd)+sc*sqrt(SX(ind(i),1:indEnd)),fliplr(muX(ind(i),1:indEnd)-sc*sqrt(SX(ind(i),1:indEnd)))];
            yp(isnan(yp)) = VBA_nanmean(yp);
            xp = [dTime,fliplr(dTime)];
            %hf(i) = fill(xp,yp,'r','parent',haf,'facecolor',col,'edgealpha',0,'facealpha',0.25);
            set(hf(i), 'Vertices', [xp', yp']);
        end
    end
    %end
    set(haf,'ygrid','on')
    axis(haf,'tight')
    set(haf,'xlim',[dTime(1),dTime(end)],'nextplot','add');
    if ~isempty(color)
        set(hp,'color',color)
        set(hf,'FaceColor',color)
    end

else
    
    hp = findobj(haf,'Tag','putBar');
    if isempty (hp)
        hp = bar(dTime:dTime+n-1,muX(ind),'facecolor',[.8 .8 .8],'parent',haf,'Tag','putBar');
        set(haf,'nextplot','add')
    else
        set(hp, ...
            'XData', dTime:dTime+n-1, ...
            'YData', muX(ind)); 
    end
    
    hf = findobj(haf,'Tag','putError');
    if isempty (hf)
        hf = errorbar(dTime:dTime+n-1,muX(ind),sc*sqrt(SX(ind)),'r.','parent',haf, 'Tag','putError');
    else
        set(hf, ...
            'XData', dTime:dTime+n-1, ...
            'YData', muX(ind), ...
            'LData', full(sc*sqrt(SX(ind))), ...
            'UData', full(sc*sqrt(SX(ind))));     
    end
    if ~isempty(color)
        set(hp,'FaceColor',color)
        set(hf,'color',color)
    end
    set(haf,'xlim',[0.2,n+0.8],'xtick',1:n,'nextplot','add');
end

% Add confidence intervals scaling control
try % only when no parent handle has been specified
    noButton;
    ud.mu = muX;
    ud.var = SX;
    ud.hf = hf;
    ud.ind = ind;
    hh(1) = uicontrol('style','edit','callback',@dox,'userdata',ud,'string','1','tooltipstring','X times posterior standard deviation...');
    set(hh(1),'units','normalized')
    pos = get(hh(1),'position');
    set(hh(1),'position',[0.4 0.02 pos(3:4)])
    hh(2) = uicontrol('style','text','string','Change error bars','units','normalized','position',[0.5 0.02 0.2 0.0476],'backgroundcolor',ones(1,3));
end


function [] = dox(e1,e2)
ud = get(e1,'userdata');
scale = str2double(get(e1,'string'));
mu = ud.mu;
SX = ud.var;
ind = ud.ind;
[n,indEnd] = size(mu);
if indEnd > 1
    for i = ind
        yp = [mu(ind(i),:)+scale.*sqrt(SX(ind(i),:)),fliplr(mu(ind(i),:)-scale.*sqrt(SX(ind(i),:)))];
        set(ud.hf(i),'ydata',yp)
    end
else
    set(ud.hf,'LData',scale.*sqrt(SX),'UData',scale.*sqrt(SX))
end
