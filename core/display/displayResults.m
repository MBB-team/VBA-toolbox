function [hres,hres2] = displayResults(posterior,out,y,x,x0,theta,phi,alpha,sigma)
% compares VB posteriors with simulated parameters and hidden states

if isempty(posterior)
    return
end

hres = figure;
pos = get(hres,'position');
set(hres,'name','Simulation results','position',[pos(1),pos(2)-pos(4),pos(3),2*pos(4)],'color',ones(1,3),'menubar','none');

% parameters
    hs = subplot(2,2,1,'parent',hres);
    if isempty(theta)
    placeHolder(hs,'no evolution parameters')
    else
    xtick = 1:out.dim.n_theta;
    set(hs,'xtick',xtick,'nextplot','add','xlim',[.2,out.dim.n_theta+.8])
    if ~out.options.OnLine
        V = VBA_getVar(posterior.SigmaTheta);
        muTheta = posterior.muTheta;
    else
        V = VBA_getVar(posterior.SigmaTheta{end});
        muTheta = posterior.muTheta(:,end);
    end
    plotUncertainTimeSeries(muTheta,V,[],hs);
    VBA_title(hs,'theta')
    plot(hs,theta,'go')
    end

hs = subplot(2,2,2,'parent',hres);
if isempty(phi)
      placeHolder(hs,'no observation parameters')  
else
    xtick = 1:out.dim.n_phi;
    set(hs,'xtick',xtick,'nextplot','add','xlim',[.2,out.dim.n_phi+.8])
    if ~out.options.OnLine
        V = VBA_getVar(posterior.SigmaPhi);
        muPhi = posterior.muPhi;
    else
        V = VBA_getVar(posterior.SigmaPhi{end});
        muPhi = posterior.muPhi(:,end);
    end
    plotUncertainTimeSeries(muPhi,V,[],hs);
    VBA_title(hs,'phi')
    plot(hs,phi,'go')
end

% hyperparameters
hs = subplot(2,2,3,'parent',hres);
n_gs = sum([out.options.sources(:).type]==0);
if n_gs == 0
    placeHolder(hs,'no observation hyperparameters')  
else
    sigmaHat = posterior.a_sigma./posterior.b_sigma;
    vs = posterior.a_sigma./(posterior.b_sigma.^2);
    lvs = log(sigmaHat+sqrt(vs)) - log(sigmaHat);
    try
        alphaHat = posterior.a_alpha(end)./posterior.b_alpha(end);
        va = posterior.a_alpha(end)./(posterior.b_alpha(end).^2);
        lva = log(alphaHat+sqrt(va)) - log(alphaHat);
        lm = log([alphaHat;sigmaHat]);
        lv = [lva;lvs];
        xtick = [1,(1:n_gs)+1];
        
        lab = {'alpha'};
        
    catch
        lm = log(sigmaHat);
        lv = lvs;
        xtick = 1:n_gs;
        lab={};
    end
    for i=1:n_gs, lab{end+1} = ['sigma' num2str(i)]; end
    set(hs,'nextplot','add')
    plotUncertainTimeSeries(lm,lv.^2,[],hs);
    plot(hs,log([alpha;sigma(:)]),'go')
    set(hs,'xtick',xtick,'xticklabel',lab)
    ylabel(hs,'log-precisions')
    VBA_title(hs,'precision hyperparameters')
end

if out.dim.n_t > 1
    dTime = [1:out.dim.n_t];
else
    dTime = [1:out.dim.p];
end
if out.dim.n > 0
    try
        % initial conditions
        hs = subplot(2,2,4,'parent',hres);
        xtick = 1:out.dim.n;
        set(hs,'xtick',xtick,'nextplot','add','xlim',[.2,out.dim.n+.8])
        V = VBA_getVar(posterior.SigmaX0);
        plotUncertainTimeSeries(posterior.muX0,V,[],hs);
        VBA_title(hs,'initial conditions')
        plot(x0,'go')
        % Hidden states
        hres2 = figure(...
            'name','Simulation results',...
            'color',[1 1 1],...
            'menubar','none');
        hs = subplot(2,2,1,'parent',hres2);
        set(hs,'nextplot','add')
        plot(hs,x','--')
        if isempty(posterior.SigmaX)
            plot(hs,posterior.muX')
        else
            V = VBA_getVar(posterior.SigmaX.current);
            plotUncertainTimeSeries(posterior.muX,V,dTime,hs);
        end
        grid(hs,'on')
        VBA_title(hs,'estimated hidden-states time series')
    end
end

try
    hres2;
catch
    hres2 = figure(...
        'name','Simulation results',...
        'color',[1 1 1],...
        'menubar','none');
end
hs = subplot(2,2,2,'parent',hres2);
set(hs,'nextplot','add')
plot(hs,y','--')
if out.dim.n_t > 1
    plotUncertainTimeSeries(out.suffStat.gx,out.suffStat.vy,dTime,hs);
else
    plotUncertainTimeSeries(out.suffStat.gx',out.suffStat.vy',dTime,hs);
end
grid(hs,'on')
VBA_title(hs,'predicted y')

if out.dim.n > 0
    miX = min([posterior.muX(:);x(:)]);
    maX = max([posterior.muX(:);x(:)]);
    hs = subplot(2,2,3,'parent',hres2);
    set(hs,'nextplot','add')
    plot(hs,posterior.muX(:),x(:),'.')
    plot(hs,[miX,maX],[miX,maX],'r')
    VBA_title(hs,'x(t) vs <x(t)>')
    grid(hs,'on')
    axis(hs,'tight')
end
miy = min([out.suffStat.gx(:);y(:)]);
may = max([out.suffStat.gx(:);y(:)]);
hs = subplot(2,2,4,'parent',hres2);
set(hs,'nextplot','add')
plot(hs,[miy,may],[miy,may],'r')
if n_gs>0
    plot(hs,out.suffStat.gx(:),y(:),'.')
else
    [stacky,stdy,gridg] = VBA_Bin2Cont(out.suffStat.gx,y);
    errorbar(gridg,stacky,stdy,'k.','parent',hs)
end

VBA_title(hs,'y(t) vs <y(t)>')
grid(hs,'on')
axis(hs,'tight')

try, VBA_getSubplots (); end

end

% display low key text when usual plot is not required
 function placeHolder(h,label)
     xx = get(h,'XLim');
     yy = get(h,'YLim');
     t=text(mean(xx),mean(yy),label,'parent',h);
     set(t, ...
         'HorizontalAlignment','center'     , ...
         'FontSize'           ,10           , ...
         'Color'              ,[.6 .6 .6]   );
    set(h,'Visible','off');
 end
