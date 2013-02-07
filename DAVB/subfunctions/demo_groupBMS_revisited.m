% Polya's urn

% close all
clear all
clc

gridn = 2.^[4,6]; % # subjects
K = 3; % # models
f = ones(K,1)./K; % true frequency of models
a0 = 1; % priors counts for each model
alpha = 0.05; % significance level
N = 256; % # Monte-Carlo simulations

hf = figure('color',[1 1 1],'name',['K=',num2str(K)]);
str = cell(length(gridn),1);
for i=1:length(gridn)
    
    str{i} = ['n=',num2str(gridn(i))];
    
    ha(i,1) = subplot(length(gridn),3,(i-1)*3+1,'parent',hf,'nextplot','add','xlim',[0,1],'ylim',[0,1]);
    title(ha(i,1),str{i})
    xlabel(ha(i,1),'Bayes Factor')
    ylabel(ha(i,1),'exceedance prob')
    plot(ha(i,1),[0,1],[1,1]./K,'r--')
    maHa = (gridn(i)/K+a0)*ones(K,1); % max-entropic marbles sample
    miHa = [gridn(i);zeros(K-1,1)]+a0; % min-entropic marbles sample
    maF = gammaln(K*a0) - K*gammaln(a0) + sum(gammaln(miHa)) - gammaln(gridn(i)+K*a0);
    miF = gammaln(K*a0) - K*gammaln(a0) + sum(gammaln(maHa)) - gammaln(gridn(i)+K*a0);
    F0 = -gridn(i)*log(K);
    mir = sigm(miF-F0);
    mar = sigm(maF-F0);
    plot(ha(i,1),[mir,mir],[0,1],'r--')
    plot(ha(i,1),[mar,mar],[0,1],'r--')
    
    ha(i,2) = subplot(length(gridn),3,(i-1)*3+2,'parent',hf,'nextplot','add','xlim',[0,1],'ylim',[0,1]);
    title(ha(i,2),str{i})
    xlabel(ha(i,2),'Bayes Factor')
    ylabel(ha(i,2),'estimated freq')
    mif = 1./K;
    maf = (gridn(i)+a0)./(gridn(i)+K*a0);
    plot(ha(i,2),[0,1],[mif,mif],'r--')
    plot(ha(i,2),[0,1],[maf,maf],'r--')
    plot(ha(i,2),[mir,mir],[0,1],'r--')
    plot(ha(i,2),[mar,mar],[0,1],'r--')
    
    ha(i,3) = subplot(length(gridn),3,(i-1)*3+3,'parent',hf,'nextplot','add','xlim',[0,1],'ylim',[0,1]);
    title(ha(i,3),str{i})
    xlabel(ha(i,3),'estimated freq')
    ylabel(ha(i,3),'exceedance prob')
    plot(ha(i,3),[0,1],[1,1]./K,'r--')
    plot(ha(i,3),[mif,mif],[0,1],'r--')
    plot(ha(i,3),[maf,maf],[0,1],'r--')
end

FPR = zeros(length(gridn),2);
r = zeros(length(gridn),N);
mep = zeros(length(gridn),N);
mEf = zeros(length(gridn),N);
for i=1:length(gridn)
    n = gridn(i);
    for j=1:N
        m = sampleFromArbitraryP(f,eye(length(f)),n)';
        [F,F0,a] = PolyaLEV(m,a0);
        ep = VBA_ExceedanceProb(a,[],'dirichlet',0);
        Ef = a./sum(a);
        r(i,j) = sigm(F-F0);
        mep(i,j) = max(ep);
        mEf(i,j) = max(Ef);
        plot(ha(i,1),r(i,j),mep(i,j),'r.')
        plot(ha(i,2),r(i,j),mEf(i,j),'r.')
        plot(ha(i,3),mEf(i,j),mep(i,j),'r.')
        drawnow
    end
    
end

col = getColors(length(gridn));
hf(2) = figure('color',[1 1 1],'name',['K=',num2str(K)]);

hb(1) = subplot(3,1,1,'parent',hf(2),'nextplot','add');
[eh,gy] = empiricalHist(r');
plot(hb(1),gy,eh,'.')
for i=1:length(gridn)
    t(i) = gy(min(find(cumsum(eh(:,i))>=1-alpha)),i);
    plot(hb(1),[t(i),t(i)],[0,max(eh(:))],'color',col(i,:))
    maHa = (gridn(i)/K+a0)*ones(K,1); % max-entropic marbles sample
    miHa = [gridn(i);zeros(K-1,1)]+a0; % min-entropic marbles sample
    maF = gammaln(K*a0) - K*gammaln(a0) + sum(gammaln(miHa)) - gammaln(gridn(i)+K*a0);
    miF = gammaln(K*a0) - K*gammaln(a0) + sum(gammaln(maHa)) - gammaln(gridn(i)+K*a0);
    F0 = -gridn(i)*log(K);
    mir = sigm(miF-F0);
    mar = sigm(maF-F0);
    plot(hb(1),[mir,mir],[0,max(eh(:))],'color',col(i,:),'linestyle','--')
    plot(hb(1),[mar,mar],[0,max(eh(:))],'color',col(i,:),'linestyle','--')
end
set(hb(1),'xlim',[0,1],'ylim',[0,max(eh(:))])
title(hb(1),'BF''s histogram')
legend(hb(1),str)

hb(2) = subplot(3,1,2,'parent',hf(2),'nextplot','add');
[eh,gy] = empiricalHist(mep');
plot(hb(2),gy,eh,'.')
for i=1:length(gridn)
    t(i) = gy(min(find(cumsum(eh(:,i))>=1-alpha)),i);
    plot(hb(2),[t(i),t(i)],[0,max(eh(:))],'color',col(i,:))
    mif = 1./K;
    plot(hb(2),[mif,mif],[0,max(eh(:))],'k--')
end
set(hb(2),'xlim',[0,1],'ylim',[0,max(eh(:))])
title(hb(2),'max(ep)''s histogram')

hb(3) = subplot(3,1,3,'parent',hf(2),'nextplot','add');
[eh,gy] = empiricalHist(mEf');
plot(hb(3),gy,eh,'.')
for i=1:length(gridn)
    t(i) = gy(min(find(cumsum(eh(:,i))>=1-alpha)),i);
    plot(hb(3),[t(i),t(i)],[0,max(eh(:))],'color',col(i,:))
    mif = (ceil(gridn(i)/K)+a0)./(gridn(i)+K*a0);
    maf = (gridn(i)+a0)./(gridn(i)+K*a0);
    plot(hb(3),[mif,mif],[0,max(eh(:))],'color',col(i,:),'linestyle','--')
    plot(hb(3),[maf,maf],[0,max(eh(:))],'color',col(i,:),'linestyle','--')
end
set(hb(3),'xlim',[0,1],'ylim',[0,max(eh(:))])
title('max(Ef)''s histogram')




getSubplots