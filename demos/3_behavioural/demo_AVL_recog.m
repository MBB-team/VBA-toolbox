% OTO: demo 'observing the observer'

clear variables
close all

% parameters of the perceptual model:
% theta(1): within-class stimuli precision
% theta(2): (static) volatility of cue-outcome association probability
% mu: class-specific 1st order moment of stimuli
% phi: response model parameters
% sigma: reaction-times precision
% flag: perceptual model (1: static, 2: dynamic, 3: volatile)
theta = [2;-1];
mu = [1,-1];
phi = [1;-1;2];
sigma = 1e1;
flag = 2;

n_t = 300; % number of trials


switch flag
    
    case 1
        
        % assume static cue-outcome association
        x2 = randn.*ones(1,n_t);
        sx2 = VBA_sigmoid(x2);
        sx2 = sx2(:)';
        
    case 2
        
        % sample the cue-outcome association according to AR(1) model
        f_fname = @f_AR;
        g_fname = @g_sigm;
        [sx2,x2,x20,eta,e] = VBA_simulate (n_t,f_fname,g_fname,[],[],[],exp(-theta(2)),Inf,[],0);
        
    case 3
        
        % First sample the association volatility according to AR(1) model
        f_fname = @f_AR;
        g_fname = @g_Id;
        [x3,x3,x30,eta,e] = VBA_simulate (n_t,f_fname,g_fname,[],[],[],exp(-theta(2)),Inf,[],0);
        % create prior variance structure for cue-outcome association
        ex3 = exp(-x3);
        for t=1:n_t
            opt.priors.iQx{t} = ex3(t);
        end
        % sample the cue-outcome association according to AR(1) model with
        % varying prior variance
        f_fname = @f_AR;
        g_fname = @g_sigm;
        [sx2,x2,x20,eta,e] = VBA_simulate (n_t,f_fname,g_fname,[],[],[],1,Inf,opt,0);
        
        
end

% Then sample outcome identity from binomial pdf:
x1 = zeros(1,n_t);
seed = 1e4*rand;
x1 = VBA_random ('Bernoulli', sx2);

% Finally, sample visual outcome from GMM ...
u = zeros(2,n_t);
for t=1:n_t
    if x1(t)==1
        u(1,t) = mu(1) + exp(-theta(1))*randn;
    else
        u(1,t) = mu(2) + exp(-theta(1))*randn;
    end
end
% ... and define observer's choices
u(2,:) = x1;
u(2,50:50:end) = 1-u(2,50:50:end); % -> categorization errors

%-----
% Now invert perceptual model using optimal VB observer algo:
inF.flag = flag;    % perceptual model
inF.mu = mu;
inF.n = 1;          % max # iterations (per trial) for the VB observer
inF.tdf = 1e-2;
inF.uu = 1;         % index of sensory signals in the vector u
inG.uc = 2;         % index of observer's choices in the vector u
f_fname = @f_AVL;
g_fname = @g_AVL;
options.inF = inF;
options.inG = inG;
if ismember(options.inF.flag,[1,2])
    X0 = [0.5;0;1e1;0];
elseif options.inF.flag == 3
    X0 = [0.5;0;1e1;0;1e1;0];
end
[RT,x,x0,eta,e] = VBA_simulate (n_t,f_fname,g_fname,theta,phi,u,Inf,Inf,options,X0);

hf = figure('color',[1 1 1]);
ha = subplot(2,2,1,'parent',hf,'nextplot','add');
plot(ha,u(1,:),'k.')
grid(ha,'on')
axis(ha,'tight')
title(ha,'simulated sensory signals')
ha = subplot(2,2,2,'parent',hf,'nextplot','add');
if ismember(options.inF.flag,[1,2])
    m = x(1:2,:);
    v = [x(1,:).*(1-x(1,:));x(3,:)];
elseif inF.flag == 3
    m = x([1:2,4],:);
    v = [x(1,:).*(1-x(1,:));x(3,:);x(5,:)];
end
plotUncertainTimeSeries(m,v,1:n_t,ha)
if inF.flag ~= 3
    plot(ha,[x1;x2]','.')
    legend({'cue identity','cue-outcome association'})
else
    plot(ha,[x1;x2;x3]','.')
    legend({'cue identity','cue-outcome association','volatility'})
end
grid(ha,'on')
axis(ha,'tight')
title(ha,'simulated VB observer')
ha = subplot(2,2,3,'parent',hf,'nextplot','add');
plot(ha,RT,'k.')
grid(ha,'on')
axis(ha,'tight')
title(ha,'simulated reaction times')
ha = subplot(2,2,4,'parent',hf,'nextplot','add');
miy = min([x2,x(2,:)]);
may = max([x2,x(2,:)]);
plot(ha,x2(:),x(2,:),'.')
plot(ha,[miy,may],[miy,may],'r')
grid(ha,'on')
axis(ha,'tight')
title(ha,'simulated vs observed cue-outcome association')
VBA_getSubplots ();


% Now OBSERVE THE OBSERVER:
y = RT + 1./sigma*randn(size(RT)); % add (motor) noise
% y = y  - min(y);
options.inF.flag = flag;
switch options.inF.flag
    case 1
        priors.muX0 = [0.5;0;1e0;0];
        priors.SigmaX0 = 1e0*eye(4);
        priors.SigmaX0(3,3) = 1e0;
        priors.muTheta = [theta(1);0];
        priors.SigmaTheta = 0.*1e0*eye(2);
        %         priors.SigmaTheta(1,1) = 0;
    case 2
        priors.muX0 = [0.5;0;1e1;0];
        priors.SigmaX0 = 1e0*eye(4);
        priors.muTheta = [theta(1);0];
        priors.SigmaTheta = 1e0*eye(2);
        priors.SigmaTheta(1,1) = 0;
    case 3
        priors.muX0 = [0.5;0;1e0;-2;1e0;0];
        priors.SigmaX0 = 1e0*eye(6);
        %         priors.SigmaX0(4,4) = 1;
        theta(2) = -32;
        priors.muTheta = [0;-2];
        priors.SigmaTheta = 1e0*eye(2);
        priors.SigmaTheta(1,1) = 0;
end
priors.muPhi = [0;0;0];
priors.SigmaPhi = 1e0*eye(length(priors.muPhi));
priors.a_alpha = Inf;
priors.b_alpha = 0;
priors.a_sigma = 1e0;
priors.b_sigma = 1e0;

% Build options and dim structures for model inversion
options.priors = priors;
dim.n_phi = length(priors.muPhi);
dim.n_theta = 2;
switch options.inF.flag
    case 1
        dim.n = 4;
    case 2
        dim.n = 4;
    case 3
        dim.n = 6;
end
[posterior,out] = VBA_NLStateSpaceModel(y,u,f_fname,g_fname,dim,options);


% Display results
displayResults(posterior,out,y,x,x0,theta,phi,Inf,sigma)
