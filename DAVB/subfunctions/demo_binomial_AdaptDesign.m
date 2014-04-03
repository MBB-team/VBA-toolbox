% demo for binomial data inversion with adaptative design
% This demo simulates a psychophysics paradigm similar to a signal
% detection task, whereby the detection probability is a sigmoidal function
% of the stimulus contrast (which is the design control variable). However,
% neither does one know the inflexion point (detection threshold) nor the
% sigmoid steepness (d prime). Thus, the design is adpated online, in the
% aim of providing the most efficient estimate of these model parameters,
% given trial-by-trial subjects' binary choice data (seen/unseen).

clear variables
close all


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
options.binomial = 1;
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
displayResults(posterior,out,y,[],[],[],phi,[],[])


% summarize results of adaptive design strategy
[handles] = displayUncertainSigmoid(posterior,out);
set(handles.ha0,'nextplot','add')
qx = g_sigm_binomial([],phi,gridu,[]);
plot(handles.ha0,gridu,qx,'k--')
VBA_ReDisplay(posterior,out)
hf = figure('color',[1 1 1]);
ha = axes('parent',hf);
plot(ha,eu,'k','marker','.');
ylabel(ha,'design efficiency')
xlabel(ha,'trials')
box(ha,'off')
set(ha,'ygrid','on')
