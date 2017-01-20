% see demo_sparsePriors

close all
clear all

% 0- evaluate the "sparsify" function
x = -4:0.01:4;
gridP = [-log(2),0,log(2)];
for i = 1:length(gridP)
    P = gridP(i);
    sx(i,:) = sparsify(x,P);
    isx(i,:) = invSparsify(x,P);
end
hf = figure('color',[1 1 1],'name','sparsify mappings');
ha = subplot(2,2,1,'parent',hf,'nextplot','add');
plot(ha,x,sx')
legend(ha,{'sparsity exponent=1/2','sparsity exponent=1','sparsity exponent=2'})

ha = subplot(2,2,2,'parent',hf,'nextplot','add');
plot(ha,x,isx')

ha = subplot(2,2,3,'parent',hf,'nextplot','add');
plot(ha,x,-0.5*isx.^2')

ps = exp(-0.5*isx.^2);
ps = ps./repmat(sum(ps,2),1,size(ps,2));
ha = subplot(2,2,4,'parent',hf,'nextplot','add');
plot(ha,x,ps)


% 1- simulate "sparse" GLM and invert
N1 = 16;
N2 = 2;
p = 16;
n = 32;
sigma = 1;

dims.n = 0;
dims.n_theta = 0;
options.inG.sparseP = 1;
options.checkGrads = 0;
options.MaxIter = 128;
options.DisplayWin = 0;
% options.priors.muPhi = 1e-2*ones(n,1);

for i=1:N1
    for j=1:N2
        
        % simulate data
        A = randn(p,n);
        X = 2*ones(n,1);
        phi1 = X;
        phi1(1:i+4) = 0;
        y1 = A*phi1 + sqrt(sigma.^-1)*randn(p,1);
        
        % L2 estimator
        options.inG.X = A;
        dims.n_phi = n;
        [p1,o1] = VBA_NLStateSpaceModel(y1,[],[],@g_GLM,dims,options);
        Phat = p1.muPhi;
        tmp = corrcoef(phi1,Phat);
        r(1,i,j) = tmp(2,1);
        e(1,i,j) = SSE(phi1,Phat);
        
        % L1 estimator
        dims.n_phi = n;
        [p2,o2] = VBA_NLStateSpaceModel(y1,[],[],@g_GLMsparse,dims,options);
        Phat = sparseTransform(p2.muPhi,options.inG.sparseP);
        tmp = corrcoef(phi1,Phat);
        r(2,i,j) = tmp(2,1);
        e(2,i,j) = SSE(phi1,Phat);
        
        % adaptive sparse-estimator
        dims.n_phi = n+1;
        [p3,o3] = VBA_NLStateSpaceModel(y1,[],[],@g_GLMsparse2,dims,options);
        Phat = sparsify(p3.muPhi(1:n),p3.muPhi(n+1));
        tmp = corrcoef(phi1,Phat);
        r(3,i,j) = tmp(2,1);
        se(i,j) = p3.muPhi(n+1);
        e(3,i,j) = SSE(phi1,Phat);
        
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



