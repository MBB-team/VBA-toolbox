% this demo reproduces the dissociation between fit accuracy and
% generalization accuracy.
% In brief, it simulates data under a GLM compriing n=16 regressors (with
% random weights). This simulated data is then splitted into train (2/3)
% and test (1/3) datasets. GLMs equipped with 0<i<32 regressors are then
% fitted on the train set, and the accuracy of their out-of-sample
% prediction is then evaluated on the test set (test-likelihood). Different
% in-sample model acuracy are then compared, as model complexity increases.
% NB: average model fit measures are obtained from N=32 Monte-Carlo
% simulations.

clear all
close all

n = 16;
b = [ones(n,1);zeros(n,1)];

C = 0.5*ones(2*n,2*n) + 0.5*eye(2*n);

N = 20; % # Monte-Carlo samples
for ii=1:N
    
    % simulate data
    X = randn(3*n,2*n)*C;
    e = randn(3*n,1);
    Y = X*b + e;
    Y_train = Y(1:2*n);
    Y_test = Y(2*n+1:3*n);
    X_train = X(1:2*n,:);
    X_test = X(2*n+1:3*n,:);
    
    for i = 1:2*n
        
        % fit model on train set
        options.DisplayWin = 0;
        options.verbose = 0;
        options.inG.X = X_train(:,1:i);
        options.priors = [];
        options.priors.SigmaPhi = 1e2*eye(i);
        dim.n = 0;
        dim.n_theta = 0;
        dim.n_phi = i;
        dim.p = 2*n;
        dim.n_t = 1;
        [p,o] = VBA_NLStateSpaceModel(Y_train,[],[],@g_GLM,dim,options);
        R2(i,ii) = o.fit.R2;
        LL(i,ii) = o.fit.LL;
        F(i,ii) = o.F;
        AIC(i,ii) = o.fit.AIC;
        BIC(i,ii) = o.fit.BIC;
        
        % derive out-of-samle prediction on test set and measure generalization error
        dim.p = n;
        dim.n_t = 1;
        options.inG.X = X_test(:,1:i);
        options.priors = p;
        options.priors.iQy{1,1} = eye(n,n);
        [muy,Vy] = VBA_getLaplace([],[],@g_GLM,dim,options);
        ge(i,ii) = sum((Y_test-muy).^2);
        
    end
end
        
hf = figure('color',[1 1 1]);
ha = axes('parent',hf,'nextplot','add');
LLg = -0.5*ge - 0.5*n*log(2*pi);
errorbar(mean(LLg'),std(LLg')./sqrt(N),'parent',ha,'marker','.');
errorbar(mean(LL'),std(LL')./sqrt(N),'parent',ha,'marker','.','color','r');
errorbar(mean(F'),std(F')./sqrt(N),'parent',ha,'marker','.','color','g');
errorbar(mean(AIC'),std(AIC')./sqrt(N),'parent',ha,'marker','.','color','y');
errorbar(mean(BIC'),std(BIC')./sqrt(N),'parent',ha,'marker','.','color','m');
plot([n,n],get(ha,'ylim'),'k--')
legend(ha,{'test Log-Likelihood','train Log-Likelihood','Free Energy','AIC','BIC'})
ylim([-500, 0])


