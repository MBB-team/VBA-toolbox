function [p] = VBA_PPM(m,v,t,disp,form)
% computes the exceedance probability of a random effect (+ plots)
% FORMAT [p] = VB_PPM(m,v,t,disp,form)
% IN:
%   - m/v: sufficient statistics of the pdf
%       -> if form='gaussian', m=E[x] and v=V[x]
%       -> if form='beta', m=counts for y=1 and v=counts for y=0, where y
%       would be the conjugate binomial variable.
%   - t: if t is a scalar, then it is a threshold (to calculate p =
%   P(x>t)) ; if it is a vector, then it is an interval (to calculate p =
%   P(x<t(1) or x>t(2)).
%   - disp: a binary switch that flags whether or not to display the
%   exceedance probability onto the posterior pdf.
%   - form: 'gaussian' or 'beta'
% OUT:
%   - p: the exceedance probability

try, disp; catch disp=0; end
try, form; catch form='gaussian'; end

switch form
    case 'gaussian'
        s = sqrt(v);
        dx = s*1e-4;
        gridx = -8*s:dx:8*s;
        gridx = gridx + m;
        f = exp(-0.5*(gridx-m).^2./v);
    case 'beta'
        gridx = 0:1e-2:1;
        logf = (m-1).*log(gridx) + (v-1).*log(1-gridx);
        f = exp(logf-max(logf));
%         f = (gridx.^(m-1)).*((1-gridx).^(v-1));
end
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
    ha = axes('parent',hf,'nextplot','add');
    if length(t) == 1
        yp = [f(indt:end),zeros(size(f(indt:end)))];
        xp = [gridx(indt:end),fliplr(gridx(indt:end))];
        fill(xp,yp,'r','parent',ha,'facecolor',0.5*[1 0 0],'edgealpha',0,'facealpha',0.25);
        plot(ha,[gridx(indt),gridx(indt)],[0,f(indt)],'r--')
        title(ha,['P(x>',num2str(t,'%4.3e'),') = ',num2str(p,'%4.3e')])
    else
        yp1 = [f(1:indt1),zeros(size(f(1:indt1)))];
        xp1 = [gridx(1:indt1),fliplr(gridx(1:indt1))];
        yp2 = [f(indt2:end),zeros(size(f(indt2:end)))];
        xp2 = [gridx(indt2:end),fliplr(gridx(indt2:end))];
        fill(xp1,yp1,'r','parent',ha,'facecolor',0.5*[1 0 0],'edgealpha',0,'facealpha',0.25);
        fill(xp2,yp2,'r','parent',ha,'facecolor',0.5*[1 0 0],'edgealpha',0,'facealpha',0.25);
        plot(ha,[gridx(indt1),gridx(indt1)],[0,f(indt1)],'r--')
        plot(ha,[gridx(indt2),gridx(indt2)],[0,f(indt2)],'r--')
        title(ha,['P(x<',num2str(t(1),'%4.3e'),' or x>',...
            num2str(t(2),'%4.3e'),') = ',num2str(p,'%4.3e')])
    end
    plot(ha,gridx,f,'k')
    grid(ha,'on')
    axis(ha,'tight')
    box(ha,'off')
    xlabel(ha,'effect of interest: x')
    ylabel(ha,'probability density function: p(x)')
    try;getSubplots;end
end








