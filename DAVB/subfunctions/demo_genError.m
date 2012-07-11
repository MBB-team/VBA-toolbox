% demo for generalization error on continuous regression

close all
clear variables

nu = 16; % number of inputs
ns = 128; % number of samples or trials
inG.dim.p = 16; % number of outputs

u = randn(nu,ns);
phi = 2e0 + randn((nu+1)*inG.dim.p,1);
sigma = 1e-1;            % precision 
g_fname = @g_matmap;

% Build simulated observations
options.inG = inG;
[y,x,x0,eta,e] = simulateNLSS(ns,[],g_fname,[],phi,u,[],sigma,options,[]);
gx = y -e;
% display time series of hidden states and observations
figure,
plot(y','ro')
hold on;
plot(gx')


% pause


%---- Invert mapping on simulated data ----%
dim.n_phi = (nu+1)*inG.dim.p; % nb of observation parameters
dim.n_theta = 0; % nb of evolution parameters
dim.n=0; % nb of hidden states
% Build priors structure

priors.muPhi = zeros(dim.n_phi,1);         % prior mean on observation params
priors.SigmaPhi = 1e0*eye(dim.n_phi); % prior covariance on observation params
priors.a_sigma = 1e0;             % Jeffrey's prior
priors.b_sigma = 1e0;             % Jeffrey's prior
% Build options structure
options.checkGrads = 0;
options.priors = priors;        % include priors in options structure
options.verbose = 1;

% leave-one-out strategy
P = zeros(dim.n_phi,ns);
my = zeros(inG.dim.p,ns);
Vy = zeros(inG.dim.p,ns);
gridalpha = 5e-2:1e-2:5e-1;
na = length(gridalpha);
e = zeros(ns,na);
for i=1:ns
    i
    yi = y(:,setdiff(1:ns,i));
    ui = u(:,setdiff(1:ns,i));
    % Call inversion routine
    [posterior,out] = VBA_NLStateSpaceModel(yi,ui,[],g_fname,dim,options);
    % derive predictive density on left-out datapoint
    opt = out.options;
    opt.priors = posterior;
    DIM = out.dim;
    DIM.n_t = 1;
    [my(:,i),V] = getLaplace(u(:,i),[],g_fname,DIM,opt);
    Vy(:,i) = diag(V);
    % store params for stability check
    P(:,i) = posterior.muPhi;
    % evaluate (conservative) error
%     displayResults(posterior,out,yi,[],[],[],phi,[],sigma)
%     pause
    Y = gx(:,i);
    MY = my(:,i);
    VY = Vy(:,i);
    ny = length(MY);
    p = zeros(ny,na);
    for j1=1:length(MY)
        for j2=1:na
            % get posterior predictive credible interval
            [x] = invnormalcdf(gridalpha(j2)/2,MY(j1),VY(j1));
            dx = x - MY(j1);
            I = [MY(j1)-dx,MY(j1)+dx];
            % check if data point is within the interval
            if ~iswithin(Y(j1),I)
                e(i,j2) = e(i,j2) +1;
            end
        end
    end
    e(i,:) = e(i,:)./ny;
end

hf = figure('color',[1 1 1]);

ha = subplot(2,2,1);
plot(ha,gx(:),my(:),'.')
hold(ha,'on')
miy = min([y(:);my(:)]);
may = max([y(:);my(:)]);
plot([miy,may],[miy,may],'r')
axis(ha,'tight')
grid(ha,'on')

title(ha,'generalization error')
xlabel(ha,'data')
ylabel(ha,'predicted data')

ha = subplot(2,2,2);
C = corrcoef(P);
imagesc(C,'parent',ha)
title(ha,'mapping stability across folds')
axis(ha,'square')
colorbar('peer',ha)

ha = subplot(2,2,3,'parent',hf);
plotUncertainTimeSeries(my,Vy,[],ha)
hold(ha,'on')
plot(ha,gx','--')

ha = subplot(2,2,4,'parent',hf);
errorbar(ha,gridalpha,mean(e,1),std(e,[],1)./sqrt(ns))
hold(ha,'on')
plot(ha,gridalpha,gridalpha,'r')
title(ha,'FPR of the predictive credible intervals')
xlabel(ha,'theoretical error rate')
ylabel(ha,'observed error rate')
axis(ha,'tight')
grid(ha,'on')

getSubplots


