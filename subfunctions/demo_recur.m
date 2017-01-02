% this script demonstrates the simulation and inversion of a k-ToM learner

level=1; % k-ToM recursive depth
HS=cat(3,[1,0;0,1],[0,1;1,0]); % game payoff matrix
inF.player=1; % player role (here: 1=seeker, 2=hider) 
inF.lev=level;
inF.game=HS;
inF.fobs=@ObsRecGen;
inF.indParev=1; % Nb para evol
inF.indParobs=2; %Nb para obs
inF.dummyPar=[1;0;0];
NtotPar=inF.indParev +inF.indParobs;
inG=inF;
inG.indlev= defIndlev(level, NtotPar);
inG.npara=NtotPar;


%% simulate sequence of k-ToM choices
phi=[0.2;0.2]; % temperature and bias
theta=2; % prior opponent's volatility
N=100; % number of trials
y2 = bernoulli(.5,N)'; % opponent's choices (here dummy binomial sampling)
y1 = NaN(1,N);
gx = NaN(1,N);
x = zeros(sizeXrec(level,NtotPar),N+1);
[x(:,1)] = f_kToM(x(:,1),theta,[],inF); %initialize hidden states
for i=1:N
    gx(i)= g_kToM(x(:,i)',phi, [],inG) ; 
    y1(i)= gx(i)>.5;
    x(:,i+1)= f_kToM(x(:,i),theta, [y2(i);y1(i)],inF);
end

% return

dim.n = sizeXrec(level,NtotPar);
dim.n_phi = 2;
dim.n_theta = 1;
options.skipf = zeros(1,N);
options.skipf(1) = 1;
options.inF = inF;
options.inG = inG;
options.binomial = 1;
options.priors.muX0 = f_kToM(zeros(dim.n,1),zeros(2,1),[],inF);
options.priors.SigmaX0 = zeros(dim.n);
f_fname = @f_kToM;
g_fname = @g_kToM;
u = [zeros(2,1),[y2;y1]];
[posterior,out] = VBA_NLStateSpaceModel(y1,u,f_fname,g_fname,dim,options);



