% Builds MCMC density of double-well system

clear variables
close all

% Choose basic settings for simulations
f_fname = @f_dbw;
g_fname = @g_Id;
u       = [];
n_t = 1e3;
deltat = 1e-2;
alpha   = 1e2;
sigma   = Inf;
theta   = [];
phi     = [];


% Build options structure for temporal integration of SDE
inF.dt = deltat;
options.inF     = inF;
options.backwardLag = 0;

% Build priors for model inversion
priors.muX0 = [-10];
priors.SigmaX0 = 0;
priors.a_alpha = 6e5;
priors.b_alpha = 1e4;
priors.a_sigma = 1e6;
priors.b_sigma = 1e1;


% Build options and dim structures for model inversion
options.priors = priors;
dim.n_theta = 0;
dim.n_phi   = 0;
dim.n       = 1;



% Build time series of hidden states and observations
N = 1e3;
X = zeros(N,n_t);
for i=1:N
    [y,x,x0,eta,e] = simulateNLSS(n_t,f_fname,g_fname,theta,phi,u,alpha,sigma,options);
    X(i,:) = x;
end

maX = max(X(:));
miX = min(X(:));
DX = maX - miX;
centres = miX:DX./100:maX;
mcmc = zeros(length(centres),n_t);
et0 = clock;
fprintf(1,'Building MCMC time-dependent density...')
fprintf(1,'%6.2f %%',0)
for t=1:n_t
    [ny,nx] = hist(X(:,t),centres);
    mcmc(:,t) = ny(:);
    if mod(100*t/n_t,10) == 0
        fprintf(1,'\b\b\b\b\b\b\b\b')
        fprintf(1,'%6.2f %%',100*t/n_t)
    end
end
% Display progress
fprintf(1,'\b\b\b\b\b\b\b\b')
fprintf(1,[' OK (took ',num2str(etime(clock,et0)),' seconds).'])
fprintf(1,'\n')

figure,imagesc(mcmc)

% display time series of hidden states and observations
% displaySimulations(y,x,eta,e)
gridx = -10:1e-3:10;
p1 = exp(-2*(gridx-2).^2);
p2 = exp(-0.125*(gridx+2).^2);
V = -log(2*p1+p2);
figure,plot(gridx,V)
% disp('--paused--')
% pause

% 
% % Call inversion routine
% % options.EvoEval = 'show_potential(posterior)';
% % options.DisplayWin = 0;
% % options.gradF = 1;
% % [posterior,out] = VBA_onlineWrapper(y,u,f_fname,g_fname,dim,options);
% [posterior,out] = VBA_NLStateSpaceModel(y,u,f_fname,g_fname,dim,options);
% 
% 
% % Display results
% displayResults(posterior,out,y,x,x0,theta,phi,alpha,sigma)
% 
% % Make predictions
% try
%     options = out.options;
%     [xs,ys,xhat,vx,yhat,vy] = comparePredictions(...
%         n_t,theta,phi,u,alpha,sigma,options,posterior,dim);
% catch
%     disp('------!!Unable to form predictions!!------')
% end




