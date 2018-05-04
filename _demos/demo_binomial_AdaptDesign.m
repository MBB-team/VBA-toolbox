function demo_binomial_AdaptDesign
% demo for binomial data inversion with adaptative design
% This demo simulates a psychophysics paradigm similar to a signal
% detection task, whereby the detection probability is a sigmoidal function
% of the stimulus contrast (which is the design control variable). However,
% neither does one know the inflexion point (detection threshold) nor the
% sigmoid steepness (d prime). Thus, the design is adpated online, in the
% aim of providing the most efficient estimate of these model parameters,
% given trial-by-trial subjects' binary choice data (seen/unseen).

p = 1e2; %  number of trials
phi = [2;-1]; % simulated parameters: [log sigmoid slope ; inflexion point]
gridu = -2:5e-2:2; % set of potential design control variables

% configure simulation and VBA inversion 
dim.n_phi = 2;
dim.n_theta = 0;
dim.n=0;
dim.n_t = 1;
dim.p = p;
g_fname = @g_sigm_binomial;
options.sources = struct('type',1 ,'out', 1); % one binomial observation;
options.priors.muPhi = [0;0];
options.priors.SigmaPhi = 1e2*eye(2);
options.DisplayWin = 0;
options.verbose = 0;
opt = options;

% prepare graphical output window
posterior = options.priors;
hf = figure('color',[1 1 1]);
ha = subplot(2,1,1,'parent',hf);
ha2 = subplot(2,1,2,'parent',hf);
set(ha,'nextplot','add')
set(ha2,'nextplot','add')
xlabel(ha,'trials')
ylabel(ha,'sigmoid parameters')
xlabel(ha2,'u: design control variable (stimulus contrast)')
ylabel(ha2,'design efficiency')

% pre-allocate trial-dependent variables
y = zeros(p,1);
u = zeros(p,1);
sx = zeros(p,1);
eu = zeros(p,1);
mu = zeros(dim.n_phi,p);
va = zeros(dim.n_phi,p);
for t=1:p
    
    % update prior for design efficiency derivation
    dim.p = 1;
    opt.priors = posterior;
    
    if t==1
        u(1) = min(gridu);
        eu(1) = VBA_designEfficiency([],g_fname,dim,opt,u(1),'parameters');
    elseif t==2
        u(2) = max(gridu);
        eu(2) = VBA_designEfficiency([],g_fname,dim,opt,u(2),'parameters');
    else
        % find most efficient control variable
        for i=1:length(gridu)
            [e(i)] = VBA_designEfficiency([],g_fname,dim,opt,gridu(i),'parameters');
        end
        ind = find(e==max(e));
        u(t) = gridu(ind(1));
        eu(t) = e(ind(1));
        % display design eficiency as a function of control variable
        cla(ha2)
        plot(ha2,gridu,e)
        plot(ha2,gridu(ind),e(ind),'go')
        drawnow
    end
    
    % sample choice according to simulated params
    sx(t) = g_sigm_binomial([],phi,u(t),[]);
    [y(t)] = sampleFromArbitraryP([sx(t),1-sx(t)]',[1,0]',1);
    
    % invert model with all inputs and choices
    dim.p = t;
    [posterior,out] = VBA_NLStateSpaceModel(y(1:t),u(1:t),[],g_fname,dim,options);
    mu(:,t) = posterior.muPhi;
    va(:,t) = diag(posterior.SigmaPhi);
    
    % display posterior credible intervals
    if t > 1
        cla(ha)
        plotUncertainTimeSeries(mu(:,1:t),sqrt(va(:,1:t)),1:t,ha,1:2);
    end
    
end


% compare final estimates with simulations
displayResults(posterior,out,y,[],[],[],phi,[],[]);


% summarize results of adaptive design strategy
[handles] = displayUncertainSigmoid(posterior,out);
set(handles.ha0,'nextplot','add')
qx = g_sigm_binomial([],phi,gridu,[]);
plot(handles.ha0,gridu,qx,'k--')
VBA_ReDisplay(posterior,out);
hf = figure('color',[1 1 1]);
ha = axes('parent',hf);
plot(ha,eu,'k','marker','.');
ylabel(ha,'design efficiency')
xlabel(ha,'trials')
box(ha,'off')
set(ha,'ygrid','on')

end

function [handles] = displayUncertainSigmoid(posterior,out,hf)
% display analysis results

try
    handles.hf0 = hf;
catch
    handles.hf0 = figure('color',ones(1,3));
end
handles.ha0 = subplot(2,1,1,'parent',handles.hf0,'nextplot','add');
title(handles.ha0,'experimental response')
xlabel(handles.ha0,'u: design control variable (motion coherency)')
ylabel(handles.ha0,'detection probability')
grid(handles.ha0,'on')

y = out.y;
u = out.u;


% check sigmoid-like response curve
ne = 8; % number of bins
mi = min(u);
ma = max(u);
dx = (ma-mi)./ne;
edges = mi:dx:ma;
ne = length(edges)-1;
gx = zeros(ne,1);
stdg = zeros(ne,1);
gridx = zeros(ne,1);
for i=1:ne
    ind = find(u>=edges(i)&u<=edges(i+1));
    gx(i) = mean(y(ind));
    stdg(i) = std(y(ind));
    gridx(i) = mean(edges([i,i+1]));
end
errorbar(gridx,gx,stdg,'k.','parent',handles.ha0)

% add model fit
g0 = 1/2;
slope = exp(posterior.muPhi(1))./4;
vslope = slope.^2.*posterior.SigmaPhi(1,1);
options = out.options;
options.priors = posterior;
[gx,vy] = VBA_getLaplace(sort(u),[],options.g_fname,out.dim,options,0,'diag');
plotUncertainTimeSeries(gx(:)',vy(:)',sort(u)',handles.ha0);
plot(handles.ha0,sort(u),gx,'b.')
yy = [0 g0 1];
xx = (yy-g0)./slope;
vyy = ((xx.*(slope+sqrt(vslope))+g0)-yy).^2;
[haf,handles.hf,handles.hp] = plotUncertainTimeSeries(yy,vyy,xx+posterior.muPhi(2),handles.ha0);
set(handles.hf,'facecolor',[1 0 0])
set(handles.hp,'color',[1 0 0])
plot(handles.ha0,posterior.muPhi(2),g0,'go')
sip = sqrt(posterior.SigmaPhi(2,2));
plot(handles.ha0,[posterior.muPhi(2)-sip,posterior.muPhi(2)+sip],[g0 g0],'g')
legend(handles.ha0,...
    {'binned responses',...
    'sigmoid estimate',...
    '1 sigmoid std',...
    'data samples',...
    'sigmoid slope estimate',...
    '1 sigmoid slope std',...
    'inflexion point estimate',...
    '1 inflexion point std'})
box(handles.ha0,'off')
set(handles.ha0,'xgrid','off')

[ny,nx] = hist(u);
handles.ha02 = subplot(2,1,2,'parent',handles.hf0);
bar(handles.ha02,nx,ny,'facecolor',[.8 .8 .8]);
xl0 = get(handles.ha0,'xlim');
set(handles.ha02,'xlim',xl0);
xlabel(handles.ha02,'u: design control variable (stimulus contrast)')
ylabel(handles.ha02,'empirical distribution of u')
title(handles.ha02,'empirical histogram')
grid(handles.ha02,'on')
box(handles.ha02,'off')
set(handles.ha02,'xgrid','off')

try getSubplots; end

end