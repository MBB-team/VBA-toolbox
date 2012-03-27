function [p] = VB_PPM(m,v,t,disp)
% computes the exceedance probability of a random effect (+ plots)
% FORMAT [p] = VB_PPM(m,v,t,disp)
% IN:
%   - m: the posterior mean of the effect
%   - v: the posterior variance of the effect
%   - t: if t is a scalar, then it is a threshold (to calculate p =
%   P(x>t)) ; if it is a vector, then it is an interval (to calculate p =
%   P(x<t(1) or x>t(2)).
%   - disp: a binary switch that flags whether or not to display the
%   exceedance probability onto the posterior pdf.
% OUT:
%   - p: the exceedance probability
%------------------------------------------------------------
% Copyright (C) 2012 Jean Daunizeau / License GNU GPL v2
%------------------------------------------------------------

try, disp; catch disp=0; end


s = sqrt(v);
dx = s*1e-4;
gridx = -8*s:dx:8*s;
gridx = gridx + m;
f = exp(-0.5*(gridx-m).^2./v);
f = f./sum(f);
if length(t) == 1
    [mp,indt] = min(abs(gridx-t));
    p = sum(f(indt:end));
else
    [mp,indt1] = min(abs(gridx-min(t)));
    [mp,indt2] = min(abs(gridx-max(t)));
    p = 1-sum(f(indt1:indt2));
end


if disp
    hf = figure('color',[1 1 1]);
    ha = axes(...
        'parent',hf,...
        'nextplot','add');
    if length(t) == 1
        yp = [f(indt:end),zeros(size(f(indt:end)))];
        xp = [gridx(indt:end),fliplr(gridx(indt:end))];
        fill(xp,yp,'r',...
            'parent',ha,...
            'facecolor',0.5*[1 0 0],...
            'edgealpha',0,...
            'facealpha',0.25);
        plot(ha,[gridx(indt),gridx(indt)],[0,f(indt)],'r--')
        title(ha,['P(x>',num2str(t,'%4.3e'),') = ',num2str(p,'%4.3e')])
    else
        yp1 = [f(1:indt1),zeros(size(f(1:indt1)))];
        xp1 = [gridx(1:indt1),fliplr(gridx(1:indt1))];
        yp2 = [f(indt2:end),zeros(size(f(indt2:end)))];
        xp2 = [gridx(indt2:end),fliplr(gridx(indt2:end))];
        fill(xp1,yp1,'r',...
            'parent',ha,...
            'facecolor',0.5*[1 0 0],...
            'edgealpha',0,...
            'facealpha',0.25);
        fill(xp2,yp2,'r',...
            'parent',ha,...
            'facecolor',0.5*[1 0 0],...
            'edgealpha',0,...
            'facealpha',0.25);
        plot(ha,[gridx(indt1),gridx(indt1)],[0,f(indt1)],'r--')
        plot(ha,[gridx(indt2),gridx(indt2)],[0,f(indt2)],'r--')
        title(ha,['P(x<',num2str(t(1),'%4.3e'),' or x>',...
            num2str(t(2),'%4.3e'),') = ',num2str(p,'%4.3e')])
    end
    plot(ha,gridx,f,'k')
    grid(ha,'on')
    axis(ha,'tight')
    box(ha,'on')
    xlabel(ha,'effect of interest: x')
    ylabel(ha,'probability density function: p(x)')
end








