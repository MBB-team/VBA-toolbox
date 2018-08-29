% this script demonstrates the simulation and inversion of a k-ToM learner

clear all
close all
clc

K = 1; % depth of k-ToM's recursive beliefs
PERS = 0; % flag for additional perseveration component
payoffTable = cat(3,[1,0;0,1],[0,1;1,0]); % game payoff matrix (here: hide-and-seek)
role = 1;  % player's role (here: 1=seeker, 2=hider)
diluteP = 0; % partial forgetting of opponent's level (only for 2-ToM and above)?
[options,dim] = prepare_kToM(K,payoffTable,role,diluteP);

%% add-in hard-wired perseveration component
if PERS
    inG.g0 = @g_kToM; % specify observation function to PERS wrapper
    inG.indP0 = 1:2; % indices of native OBS params
    inG.indbeta = 3; % index of PERS param
    inG.in0 = options.inG; % optional inputs to native OBS function
    options.inG = inG;
    g_fname = @g_wrap_perseveration;
else
    g_fname = @g_kToM;
end
    
    
%% simulate sequence of k-ToM choices
phi = [-2;0;1]; % temperature, bias and perseveration
theta = -2; % k-ToM's prior volatility
if diluteP==1
    theta = [theta;-2];
end
N = 50; % number of trials
o = VBA_random ('Bernoulli', 0.85, 1, N); % opponent's choices (here dummy binomial sampling)
tic
a = NaN(1,N); % agent's choices
gx = NaN(1,N);
x = zeros(dim.n,N+1);
[x(:,1)] = f_kToM(options.priors.muX0,theta,[],options.inF); %initialize hidden states
for i=1:N
    if i==1;
        ug = [NaN;NaN];
    else
        ug = [o(i-1);a(i-1)];
    end
    gx(i) = feval(g_fname,x(:,i),phi,ug,options.inG) ; 
    a(i) = gx(i)>.5;
    r(i) = payoffTable(2-a(i),2-o(i),role);
    x(:,i+1) = f_kToM(x(:,i),theta,[o(i);a(i)],options.inF);
end
toc
figure,
subplot(2,1,1),plot(x')
subplot(2,1,2),plot(gx)
hold on, plot([1,N],[0.5,0.5],'color',0.2*[1 1 1])
ic = find(r==1);
plot(ic,a(ic),'g.')
ii = find(r==0);
plot(ii,a(ii),'r.')
title(['perf=',num2str(100*length(ic)./N,3),'%'])
if PERS
    hf = unwrapKTOM(x,options.inG.in0);
else
    hf = unwrapKTOM(x,options.inG);
end
set(hf,'position',[1290 266 560 857])


%% invert k-ToM model
options.skipf = zeros(1,N);
options.skipf(1) = 1;
options.sources = struct('type',1 ,'out', 1); % one binomial observation;
options.DisplayWin = 1;
options.priors.SigmaTheta = 1e2*eye(dim.n_theta); % relax evol param
options.priors.SigmaX0 = 0*eye(dim.n);
f_fname = @f_kToM;
g_fname = @g_kToM;
u = [zeros(2,1),[o;a]];
y = a(1:50);
[posterior,out] = VBA_NLStateSpaceModel(y,u,f_fname,g_fname,dim,options);

%% Display results
displayResults(posterior,out,y,x(:,1:50),options.priors.muX0,theta,phi,[],[])

if PERS
    hf = unwrapKTOM(posterior.muX,options.inG.in0);
else
    hf = unwrapKTOM(posterior.muX,options.inG);
end
set(hf,'name',[get(hf,'name'),' : VBA inversion'])


