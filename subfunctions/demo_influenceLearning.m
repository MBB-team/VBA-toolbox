% This script simulates and inverts an "influence learner"
% This model is based upon [Hampton et al., 2008].
% Let p be the agent's prediction of the other's next move, ie the
% probability that the other will pick the first alternative option. The
% "influence learning" rule can be written as follows:
% p <- p + eta*(o-p) - 2*lambda*p*(1-p)*(a-0.5*(1-beta*invsig(p)))
% where o is the other's last move, a  is the agent's last move, eta is the
% weight of the agent's prediction error, lambda is the weight of other's
% prediction error and beta is the other's temperature.

% simulation parameters
payoffTable = cat(3,[1,0;0,1],[0,1;1,0]); % game payoff matrix (here: hide-and-seek)
role = 1;  % player's role (here: 1=seeker, 2=hider) 
options.inF = {payoffTable,role};
options.inG = {payoffTable,role};
dim.n = 1;
dim.n_theta = 3;
dim.n_phi = 2;



%% simulate sequence of k-ToM choices
x0 = [0];
theta = [log(1);log(1);0];
phi=[0;0]; % temperature and bias
N=50; % number of trials
y2 = bernoulli(.5,N)'; % opponent's choices (here dummy binomial sampling)
y1 = NaN(1,N);
gx = NaN(1,N);
x = zeros(dim.n,N+1);
x(:,1) = x0; %initialize hidden states
for i=1:N
    gx(i)= g_Hampton(x(:,i)',phi, [],options.inG) ; 
    y1(i)= gx(i)>.5;
    x(:,i+1)= f_Hampton(x(:,i),theta, [y2(i);y1(i)],options.inF);
end

figure,plot(x')


options.skipf = zeros(1,N);
options.skipf(1) = 1;
options.binomial = 1;
options.SigmaTheta = 1e0*eye(dim.n_theta);
f_fname = @f_Hampton;
g_fname = @g_Hampton;
u = [zeros(2,1),[y2;y1]];
[posterior,out] = VBA_NLStateSpaceModel(y1,u,f_fname,g_fname,dim,options);



