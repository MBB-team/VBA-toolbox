function demo_sparsify()
% see demo_sparsePriors

% 0- evaluate the "sparsify" function
x = -4:0.01:4;
gridP = [-log(2),0,log(2)];
for i = 1:length(gridP)
    P = gridP(i);
    sx(i,:) = VBA_sparsifyPrior (x, P);
end
hf = figure('color',[1 1 1],'name','sparsify mappings');
ha = subplot(1,2,1,'parent',hf);
plot(ha,x,sx')
xlabel('original parameter')
ylabel('transformed parameter')
legend(ha,{'sparsity exponent=1/2','sparsity exponent=1','sparsity exponent=2'},'Location','northwest')



ha = subplot(1,2,2,'parent',hf);
for i = 1:length(gridP)
    ps = interp1(sx(i,:),normpdf(x,0,1),x);
    ps = ps / nansum(ps);
    plot(ha,x,ps)
    hold on
end
hold off
xlabel('transformed parameter')
ylabel('prior density')


% 1- simulate "sparse" GLM and invert
N1 = 16; % nb of 0 regressors
N2 = 30; % nb of MC simuls
p = 16; % nb of observations
n = 32; % nb of regressors
sigma = .1;

dims.n = 0;
dims.n_theta = 0;
options.DisplayWin = 0;
options.verbose = 0;

for i=1:N1
    for j=1:N2
        
        % simulate data
        X = randn(p,n);
        phi1 = ones(n,1);
        phi1(1:n-N1+i) = 0;
        phi1(isnan(phi1)) = 0;
        y1 = X*phi1 + sqrt(sigma)*randn(p,1);
        options.inG.X = X;

        % L2 estimator
        dims.n_phi = n;
        [p1,o1] = VBA_NLStateSpaceModel(y1,[],[],@g_GLM,dims,options);
        Phat = p1.muPhi;
        tmp = corrcoef(phi1,Phat);
        r(1,i,j) = tmp(2,1);
        e(1,i,j) = SSE(y1,o1.suffStat.gx);
        
        % L1 estimator
        dims.n_phi = n;
        [p2,o2] = VBA_NLStateSpaceModel(y1,[],[],@g_GLMsparse,dims,options);
        Phat = VBA_sparsifyPrior (p2.muPhi);
        tmp = corrcoef(phi1,Phat);
        r(2,i,j) = tmp(2,1);
        e(2,i,j) = SSE(y1,o2.suffStat.gx);
        
        % adaptive sparse-estimator
        dims.n_phi = n+1;
        options3 = options;
        options3.priors.muPhi = [zeros(n,1); log(2)];
        options3.priors.SigmaPhi = diag([ones(n,1); .2]);
        [p3,o3] = VBA_NLStateSpaceModel(y1,[],[],@g_GLMsparse2,dims,options3);
        Phat = VBA_sparsifyPrior (p3.muPhi(1:end-1), p3.muPhi(end));
        tmp = corrcoef(phi1,Phat);
        r(3,i,j) = tmp(2,1);
        se(i,j) = exp(p3.muPhi(n+1));
        e(3,i,j) = SSE(y1,o3.suffStat.gx);
        
    end
end

hf = figure('color',[1 1 1],'name','estimation accuracy');
ha = subplot(2,2,1,'parent',hf,'nextplot','add');
mr = mean(r,3);
sr = std(r,[],3);
errorbar(ha,mr',sr')
legend(ha,{'L2','L1','adaptive sparse'})
xlabel(ha,'sparsity')
ylabel(ha,'estimation accuracy')
ha = subplot(2,2,2,'parent',hf,'nextplot','add');
mse = mean(se,2);
sse = std(se,[],2);
errorbar(ha,mse,sse)
xlabel(ha,'sparsity')
ylabel(ha,'estimated sparsity exponent')
ha = subplot(2,2,3,'parent',hf,'nextplot','add');
me = mean(e,3);
se = std(e,[],3);
errorbar(ha,me',se')
xlabel(ha,'sparsity')
ylabel(ha,'SSE')

end

function sse = SSE(x,y)
% sum-of-squared distance between x and y
sse = sum(abs(vec(x)-vec(y)));
end
