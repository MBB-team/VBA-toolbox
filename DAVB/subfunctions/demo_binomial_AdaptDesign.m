% demo for binomial data inversion with adaptative design

clear variables
close all


p = 2e2; %  number of trials
phi = [2.5;-0.25]; % simulated parameters: [log sigmoid slope ; inflexion point]
gridu = -1:2e-2:1; % set of potential design control variables
optimDesign = 1; % if 1: further optimize design locally

% get binomial sampling probabilities
qx = g_sigm_binomial([],phi,gridu,[]);
hf0 = figure;
ha0 = subplot(2,1,1,'parent',hf0);
plot(ha0,gridu,qx,'k--')
xlabel(ha0,'u: design control variable (coherency)')
ylabel(ha0,'proba of choosing ''outward'' (sigmoid)')
legend(ha0,{'sampling proba: p(y|x)'})
% pause

% sample binomial data
addpath('../sampling')
addpath('../optimDesign')

y = zeros(p,1);
u = zeros(p,1);
sx = zeros(p,1);

seed = 1e4*rand;

dim.n_phi = 2;                  % nb of observation parameters
dim.n_theta = 0;                % nb of evolution parameters
dim.n=0;                        % nb of hidden states
dim.n_t = 1;
dim.p = p;

g_fname = @g_sigm_binomial;     % observation function
f_fname = [];                   % evolution function (learning)

% Call inversion routine
options.binomial = 1;
options.gradF = 0;
options.priors.muPhi = [2;0];
options.priors.SigmaPhi = 1e8*eye(2);
options.DisplayWin = 0;
options.TolFun = 1e-8;
options.GnTolFun = 1e-8;
options.Laplace = 0;
options.verbose = 0;
opt = options;

mu = zeros(dim.n_phi,p);
va = zeros(dim.n_phi,p);
posterior = options.priors;
hf = figure;
ha = subplot(2,1,1,'parent',hf);
ha2 = subplot(2,1,2,'parent',hf);
set(ha,'nextplot','add')
set(ha2,'nextplot','add')
xlabel(ha,'trials')
ylabel(ha,'sigmoid parameters')
xlabel(ha2,'u: design control variable (coherency)')
ylabel(ha2,'design efficiency')
for t=1:p
    
    % first optimize experimental design
    dim.p = 1;
    opt.priors = posterior;
    for i=1:length(gridu)
        [e(i)] = designEfficiency(f_fname,g_fname,dim,opt,gridu(i),'parameters');
    end
    ind = find(e==max(e));
    ustar = gridu(ind);
    if optimDesign
        hil = inline(...
            'designEfficiency([],args{1},args{2},args{3},x,''parameters'')',...
            'x','args');
        OPT.args = {{g_fname,dim,opt}};
        OPT.minimize = 0;
        OPT.GnMaxIter = length(gridu);
        if t>1
            u0 = u(t-1);
        else
            u0 = 0;
        end
        [ou,curv,out] = optimCost(hil,u0,OPT);
        [oe] = designEfficiency([],g_fname,dim,opt,ou,'parameters');
        u(t) = ou;
    else
        u(t) = ustar(1);
    end
    
    
    cla(ha2)
    plot(ha2,gridu,e)
    plot(ha2,gridu(ind),e(ind),'go')
    if optimDesign
        plot(ha2,ou,oe,'r+')
    end
    drawnow
%     pause(.01)
    
    % then sample choice according to simulated params
    sx(t) = g_sigm_binomial([],phi,u(t),[]);
    try
        [y(t),seed] = binomial_sample(1,sx(t),seed);
    catch
        [y(t)] = sampleFromArbitraryP([sx(t),1-sx(t)],[1,0],1);
    end
    
    % finally, invert model with all inputs and choices
    dim.p = t;
    [posterior,out] = VBA_NLStateSpaceModel(...
        y(1:t),u(1:t),[],g_fname,dim,options);
    mu(:,t) = posterior.muPhi;
    va(:,t) = diag(posterior.SigmaPhi);
    if t > 1
        cla(ha)
        plotUncertainTimeSeries(mu(:,1:t),sqrt(va(:,1:t)),1:t,ha,1:2);
%         set(ha,'xtick',1:t)
    end
    
    
    
end


% qxhat = g_sigm_binomial([],posterior.muPhi,gridu,[]);
% qx2 = g_sigm_binomial([],posterior.muPhi,u,[]);
% plot(ha0,gridu,qxhat,'k--')
% plot(ha0,u,qx2,'b.')
% legend({'sampling proba: p(y|x)','sigmoid estimate','sampled points'})


%---- Display results ----%
displayResults(posterior,out,y,[],[],[],phi,[],[])


% graphical output
set(ha0,'nextplot','add')
g0 = 1/2;
slope = exp(posterior.muPhi(1))./4;
vslope = slope.^2.*posterior.SigmaPhi(1,1);
options.priors = posterior;
dim.p = length(gridu);
[gx,vy] = getLaplace(gridu(:),[],g_fname,dim,options);
gxhat = g_sigm_binomial([],posterior.muPhi,sort(u),[]);
vy = diag(vy);
plotUncertainTimeSeries(gx(:)',vy(:)',gridu(:)',ha0);
plot(ha0,sort(u),gxhat,'b.')
plot(ha0,posterior.muPhi(2),g0,'ro')
yy = [0 g0 1];
xx = (yy-g0)./slope;
vyy = ((xx.*(slope+sqrt(vslope))+g0)-yy).^2;
[haf,hf,hp] = plotUncertainTimeSeries(yy,vyy,xx+posterior.muPhi(2),ha0);
set(hf,'facecolor',[1 0 0])
set(hp,'color',[1 0 0])
legend(ha0,...
    {'simulated response',...
    'sigmoid estimate',...
    '90% confidence interval',...
    'data samples',...
    'inflexion point',...
    'sigmoid slope'})

[ny,nx] = hist(u);
ha02 = subplot(2,1,2,'parent',hf0);
bar(ha02,nx,ny);
xl0 = get(ha0,'xlim');
set(ha02,'xlim',xl0);
xlabel(ha02,'u: design control variable (coherency)')
ylabel(ha02,'empirical distribution of u')

