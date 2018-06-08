function [p,f,gridx] = VBA_PPM(m,v,t,form,disp,ha)
% computes (and displays) P(x>t|m,v) or 1-P(t1<x<t2|m,v)
% FORMAT [p,f,gridx] = VBA_PPM(m,v,t,disp,form)
% IN:
%   - m/v: sufficient statistics of the pdf
%       -> if form='gaussian', m=E[x] and v=V[x]
%       -> if form='beta', m=counts for y=1 and v=counts for y=0, where y
%       would be the conjugate binomial variable.
%   - t: if t is a scalar, then it is a threshold (to calculate p =
%   P(x>t)) ; if it is a vector, then it is an interval (to calculate p =
%   P(x<t(1) or x>t(2)).
%   - form: 'gaussian' or 'beta'
%   - disp: a binary switch that flags whether or not to display the pdf.
%   - ha: axis' handle for display (relevant only if disp=1)
% OUT:
%   - p: P(x>t|m,v) or 1-P(t1<x<t2|m,v)
%   - f: pdf (either Gaussian or Beta) evaluated on a grid
%   - gridx: grid over which the pdf is evaluated

try, form; catch form='gaussian'; end
try, disp; catch disp=0; end

% derive probability density function
switch form
    case 'gaussian'
        s = sqrt(v);
        dx = s*1e-4;
        gridx = -8*s:dx:8*s;
        gridx = gridx + m;
        f = exp(-0.5*(gridx-m).^2./v);
    case 'beta'
        gridx = 0:1e-4:1;
        logf = (m-1).*log(gridx) + (v-1).*log(1-gridx);
        f = exp(logf-max(logf));
end
f = f./sum(f);

% derive P(x>t|m,v) or 1-P(t1<x<t2|m,v)
if length(t) == 1
    [mp,indt] = min(abs(gridx-t));
    p = sum(f(indt:end));
else
    [mp,indt1] = min(abs(gridx-min(t)));
    [mp,indt2] = min(abs(gridx-max(t)));
    p = 1-sum(f(indt1:indt2));
end

% display results
if disp
    try
        set(ha,'nextplot','add')
    catch
        hf = figure('color',[1 1 1]);
        ha = axes('parent',hf,'nextplot','add');
    end
    if length(t) == 1
        yp = [f(indt:end),zeros(size(f(indt:end)))];
        xp = [gridx(indt:end),fliplr(gridx(indt:end))];
        fill(xp,yp,'r','parent',ha,'facecolor',0.8*[1 0.75 0.75],'edgecolor','none');
        plot(ha,[gridx(indt),gridx(indt)],[0,f(indt)],'r--')
        VBA_title(ha,['P(x>',num2str(t,3),') = ',num2str(p,3)])
    else
        yp1 = [f(1:indt1),zeros(size(f(1:indt1)))];
        xp1 = [gridx(1:indt1),fliplr(gridx(1:indt1))];
        yp2 = [f(indt2:end),zeros(size(f(indt2:end)))];
        xp2 = [gridx(indt2:end),fliplr(gridx(indt2:end))];
        fill(xp1,yp1,'r','parent',ha,'facecolor',0.8*[1 0.75 0.75],'edgecolor','none');
        fill(xp2,yp2,'r','parent',ha,'facecolor',0.8*[1 0.75 0.75],'edgecolor','none');
        plot(ha,[gridx(indt1),gridx(indt1)],[0,f(indt1)],'r--')
        plot(ha,[gridx(indt2),gridx(indt2)],[0,f(indt2)],'r--')
        VBA_title(ha,['P(x<',num2str(t(1),3),' or x>',num2str(t(2),3),') = ',num2str(p,3)])
    end
    plot(ha,gridx,f,'k')
    box(ha,'off')
    xlabel(ha,'effect of interest: x')
    ylabel(ha,'probability density function: p(x)')
    try;VBA_getSubplots; end
end








