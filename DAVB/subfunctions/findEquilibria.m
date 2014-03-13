function [eq,Y,out,ha,ha2] = findEquilibria(n_t,f_fname,g_fname,theta,phi,u,alpha,sigma,options,dim,N,verbose)
% detects the number of fixed points in EGT replicator dynamics
% function [eq,X,out,ha,ha2] = findEquilibria(n_t,f_fname,g_fname,theta,phi,u,alpha,sigma,options,dim,N)
% IN:
%   - n_t ... dim: the usual entries to simulateNLSS.m
%   - N: the number of Monte-Carlo repetitions {32}
% OUT:
%   - eq: the dim.nxK matrix of ESS, where K is the number of detected
%   fixed points of the EGT replicator dynamics
%   - X: the Monte-Carlo repetitions of last time samples
%   - out: the output structure of the classifier. This is useful for
%   tracking, e.g., the size of basin of attraction of the fixed points (in
%   terms of their relative frequency across Monte-Carlo repetitions).
%   - ha/ha2: handles of axes fr display

try; N; catch; N= 32; end
try; verbose; catch; verbose = 1; end

if verbose
    hf = figure('color',[1 1 1],'name','EGT fixed points');
    ha = subplot(2,1,1,'parent',hf,'nextplot','add');
    xlabel(ha,'evolutionary time')
    ylabel(ha,'traits'' frequencies')
    title(ha,'EGT replicator dynamics')
else
    ha = [];
    ha2 = [];
end

% 1- Obtain MCMC distribution of fixed points
Y = zeros(dim.n,N);
fprintf(1,'Looking for EGT equilibria...')
fprintf(1,'%6.2f %%',0)
et0 = clock;
options.verbose = 0;
flag = 1; % flag = 1: no evolutionary force
for i=1:N
    x0 = randn(dim.n,1);
    x0 = x0 - mean(x0);
    [y] = simulateNLSS(n_t,f_fname,g_fname,theta,phi,u,alpha,sigma,options,x0);
%     [xf] = findSS(f_fname,x0,theta,u,options.inF);
%     yf = feval(g_fname,xf,phi,u,options.inG);
    if verbose
        plot(ha,y','--','linewidth',0.5)
%         plot(ha,size(y,2),yf,'o')
        drawnow
    end
    dy = repmat(mean(y,2),1,size(y,2))-y;
    if max(abs(dy(:)))>5e-2
        flag = 0;
    end
    Y(:,i) = y(:,end);
    fprintf(1,repmat('\b',1,8))
    fprintf(1,'%6.2f %%',floor(100*i/N))
end
fprintf(1,repmat('\b',1,8))
fprintf(1,[' OK (took ',num2str(etime(clock,et0)),' seconds).'])
fprintf(1,'\n')

if flag
    eq = ones(dim.n,1)./dim.n;
    out.lambda = 1;
    out.varLambda = 0;
    out.gamma = 0;
    K = 1;
    str2 = ' (no evolutionary pressure)';
else
    Ym = Y-repmat(mean(Y,2),1,N);
    % approximate on a grid to help classifier
%     Ym = approxOnGrid(Ym,0:5e-2:1);
    if isequal(Ym,zeros(size(Ym)))
        K = 1;
        out.lambda = 1;
        out.varLambda = 0;
        out.gamma = 0;
        eq = mean(Y,2);
    else
        [u,ss,Kpca,v]=PCA_MoG(Ym,0.95,0);
%         [xi,Mu,F,out,K] = VBEM_GM(u',N);
        options.verbose = 0;
        options.minSumZ = 1e-2;
        options.init = 'hierarchical';
%         options.priors.b_gamma = 1e-2*ones(2,1);
        [posterior,out] = VBA_MoG(u',2,options);
        Mu = posterior.muEta;
        K = size(posterior.z,1);
        % [handles] = plotResults(Ym,xi,Mu,F,out,K);
        eq = v(:,1:Kpca)*diag(ss(1:Kpca))*Mu + repmat(mean(Y,2),1,K);
    end
    str2 = [];
end

disp(['  --> Found ',num2str(K),' ESS.'])

if verbose
    ha2 = subplot(2,1,2,'parent',hf,'nextplot','add');
    multipie(eq+eps,1:K,ones(1,K),ha2);
    box(ha2,'off')
    set(ha2,'xtick',[1:K],'ylim',[0.4,2],'ytick',[],'xlim',[0 K+1])
    f = posterior.d./sum(posterior.d);
    for i=1:K
        str = num2str(f(i),'% 10.2f');
        text(i-0.25,1.8,str,'parent',ha2);
    end
    title(ha2,['type and frequency of EGT steady states',str2])
end

