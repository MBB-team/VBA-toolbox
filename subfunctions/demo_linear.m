function [posterior,out,betas] = demo_linear(n)
% Demo for fitting a linear model

if nargin == 0
    n=5;
end
%%
n_t = n*100 ; % numbre of observations
n_u = n; % number of regressors
% sigma = .01; % observation noise
u = randn(n_u,n_t); % random dependent variables
betas = randi(10,n_u,1) ; % weights


%%
function [gx,dgdx,dgdp] = g_linear(~,P,ut,~)
    gx = P'*ut ;
    dgdx=[];
    dgdp=ut;
end

% 
priors.muPhi = zeros(n_u,1);
priors.SigmaPhi = eye(n_u,n_u);
options.priors = priors;
options.verbose=0;
options.DisplayWin=0;

%% simulate
[y,x,x0,eta,e] = simulateNLSS(n_t,[],@g_linear,[],betas,u,1,5,options,[]);


%% invert model
dim.n_phi = n_u;                  % nb of observation parameters
dim.n_theta = 0;                % nb of evolution parameters
dim.n=0;                        % nb of hidden states
dim.u=n_u;                        % nb of inputs
dim.p=1;
[posterior,out] = VBA_NLStateSpaceModel(y,u,[],@g_linear,dim,options);

end

