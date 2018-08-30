% demo for binary data classification
% This script first simulates binary data (y) under a simple logistic-like
% model, i.e.: p(y=1) = sigmoid(X*b), where X is a feature matrix, and b
% are arbitrary weights. VBA's classifier is then applied to the simulated
% data, and performance is evaluated. In particular, cross-validation
% classification accuracy is compared to bayesian model evidence, when
% using different priors.

clear
close all

% simulate binary data
n = 32; % data sample size
p = 10; % number of features
X = randn(n,p); % feature matrix
b = 1+randn(p,1); % feature weights
e = randn(n,1); % additional noise
y = +(VBA_sigmoid(X*b+e)>0.5);

% classify data using default set-up
k = n; % number of folds (k=n: leave-one-out cross-validation)
verbose = 1; % verbose mode
options = [];
sparse = 0; % sparse mode
[all] = VBA_classification(X,y,k,verbose,options,sparse);
% evaluate estimated classification weights
displayResults(all.posterior,all.out,y,[],[],[],b,[],[])

return

% now change priors and compare class. accuracy VS bayes. model evidence
v = 10.^[-4:4];
verbose = 0;
Nmc = 32; % number of Monte-Carlo simulations
n = 16; % data sample size
p = 4; % number of features
k = n; % number of folds
for i=1:Nmc
    i
    % simulate binary data
    X = randn(n,p); % feature matrix
    b = ones(p,1); % feature weights
    e = randn(n,1); % additional noise
    y = +(VBA_sigmoid(X*b+e)>0.5);
    % apply classifier
    for j=1:length(v)
        options.priors.SigmaPhi = v(j).*eye(p);
        [all] = VBA_classification(X,y,k,verbose,options,sparse);
        F(i,j) = all.out.F;
        pa(i,j) = all.stat.bpa;
    end
end
hf = figure('color',[1 1 1],'name','classification accuracy VS model evidence');
ha = subplot(1,2,1);
errorbar(v,mean(pa,1),std(pa,[],1)./sqrt(Nmc),'parent',ha)
set(ha,'xscale','log','xlim',[10^-5,10^5],'box','off')
xlabel(ha,'priors variance')
ylabel(ha,'classification accuracy (cross-validation)')
ha = subplot(1,2,2);
errorbar(v,mean(F,1),std(F,[],1)./sqrt(Nmc),'parent',ha)
set(ha,'xscale','log','xlim',[10^-5,10^5],'box','off')
xlabel(ha,'priors variance')
ylabel(ha,'log model evidence')

