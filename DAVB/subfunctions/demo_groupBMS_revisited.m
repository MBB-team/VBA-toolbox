% Polya's urn

% close all
clear all
clc

gridn = 2.^[4,6]; % # subjects
K = 4; % # models
f = ones(K,1)./K; % true frequency of models
a0 = 1; % priors counts for each model
alpha = 0.05; % significance level
N = 2^10; % # Monte-Carlo simulations

hf = figure('color',[1 1 1],'name',['K=',num2str(K)]);
str = cell(length(gridn),1);
for i=1:length(gridn)
    
    str{i} = ['n=',num2str(gridn(i))];
    
    ha(i,1) = subplot(length(gridn),3,(i-1)*3+1,'parent',hf,'nextplot','add','xlim',[0,1],'ylim',[0,1]);
    title(ha(i,1),str{i})
    xlabel(ha(i,1),'Bayes Factor')
    ylabel(ha(i,1),'exceedance prob')
    plot(ha(i,1),[0,1],[1,1]./K,'b--')
    maHa = (gridn(i)/K+a0)*ones(K,1); % max-entropic marbles sample
    miHa = [gridn(i);zeros(K-1,1)]+a0; % min-entropic marbles sample
    maF = gammaln(K*a0) - K*gammaln(a0) + sum(gammaln(miHa)) - gammaln(gridn(i)+K*a0);
    miF = gammaln(K*a0) - K*gammaln(a0) + sum(gammaln(maHa)) - gammaln(gridn(i)+K*a0);
    F0 = -gridn(i)*log(K);
    mir = sigm(miF-F0);
    mar = sigm(maF-F0);
    plot(ha(i,1),[mir,mir],[0,1],'b--')
    plot(ha(i,1),[mar,mar],[0,1],'b--')
    
    ha(i,2) = subplot(length(gridn),3,(i-1)*3+2,'parent',hf,'nextplot','add','xlim',[0,1],'ylim',[0,1]);
    title(ha(i,2),str{i})
    xlabel(ha(i,2),'Bayes Factor')
    ylabel(ha(i,2),'estimated freq')
    mif = 1./K;
    maf = (gridn(i)+a0)./(gridn(i)+K*a0);
    plot(ha(i,2),[0,1],[mif,mif],'b--')
    plot(ha(i,2),[0,1],[maf,maf],'b--')
    plot(ha(i,2),[mir,mir],[0,1],'b--')
    plot(ha(i,2),[mar,mar],[0,1],'b--')
    
    ha(i,3) = subplot(length(gridn),3,(i-1)*3+3,'parent',hf,'nextplot','add','xlim',[0,1],'ylim',[0,1]);
    title(ha(i,3),str{i})
    xlabel(ha(i,3),'estimated freq')
    ylabel(ha(i,3),'exceedance prob')
    plot(ha(i,3),[0,1],[1,1]./K,'b--')
    plot(ha(i,3),[mif,mif],[0,1],'b--')
    plot(ha(i,3),[maf,maf],[0,1],'b--')
end

FPR = zeros(length(gridn),2);
r = zeros(length(gridn),N,2);
mep = zeros(length(gridn),N,2);
mEf = zeros(length(gridn),N,2);
for i=1:length(gridn)
    n = gridn(i);
    fprintf(1,['n= ',num2str(gridn(i)),' ...'])
    fprintf(1,'%6.2f %%',0)
    for j=1:N
        % invert model under the null
        m = sampleFromArbitraryP(f,eye(length(f)),n)';
        [F,F0,a] = PolyaLEV(m,a0);
        ep = VBA_ExceedanceProb(a,[],'dirichlet',0);
        Ef = a./sum(a);
        r(i,j,1) = sigm(F-F0);
        mep(i,j,1) = max(ep);
        mEf(i,j,1) = max(Ef);
        % invert model under the alternative
        f2 = VBA_sample('dirichlet',struct('d',a0*ones(K,1)),1,0);
        m = sampleFromArbitraryP(f2,eye(length(f)),n)';
        [F,F0,a] = PolyaLEV(m,a0);
        ep = VBA_ExceedanceProb(a,[],'dirichlet',0);
        Ef = a./sum(a);
        r(i,j,2) = sigm(F-F0);
        mep(i,j,2) = max(ep);
        mEf(i,j,2) = max(Ef);
        % display statistics
        plot(ha(i,1),r(i,j,1),mep(i,j,1),'r.')
        plot(ha(i,2),r(i,j,1),mEf(i,j,1),'r.')
        plot(ha(i,3),mEf(i,j,1),mep(i,j,1),'r.')
        plot(ha(i,1),r(i,j,2),mep(i,j,1),'g.')
        plot(ha(i,2),r(i,j,2),mEf(i,j,1),'g.')
        plot(ha(i,3),mEf(i,j,2),mep(i,j,1),'g.')
        drawnow
        if isequal(mod(j,ceil(N/200)),0)
            fprintf(1,repmat('\b',1,8))
            fprintf(1,'%6.2f %%',100*j/N)
        end
    end
    fprintf(1,repmat('\b',1,8))
    fprintf(' OK.')
    fprintf('\n')
    
end

% col = getColors(length(gridn));
hf(2) = figure('color',[1 1 1],'name',['K=',num2str(K)]);
s = 1/4;

s0 = {'H0','H1'};

col = ['r','g'];
may = zeros(3,1);
for ij=1:2
    
    hb(1) = subplot(3,1,1,'parent',hf(2),'nextplot','add');
    [eh,gy] = empiricalHist(r(:,:,ij)',s);
    plot(hb(1),gy(:,1),eh(:,1),'color',col(ij),'linestyle','-')
    plot(hb(1),gy(:,2),eh(:,2),'color',col(ij),'linestyle','.')
    %     for i=1:length(gridn)
    %         t(i) = gy(min(find(cumsum(eh(:,i))>=1-alpha)),i);
    %         plot(hb(1),[t(i),t(i)],[0,max(eh(:))],'color',col(ij))
    %         maHa = (gridn(i)/K+a0)*ones(K,1); % max-entropic marbles sample
    %         miHa = [gridn(i);zeros(K-1,1)]+a0; % min-entropic marbles sample
    %         maF = gammaln(K*a0) - K*gammaln(a0) + sum(gammaln(miHa)) - gammaln(gridn(i)+K*a0);
    %         miF = gammaln(K*a0) - K*gammaln(a0) + sum(gammaln(maHa)) - gammaln(gridn(i)+K*a0);
    %         F0 = -gridn(i)*log(K);
    %         mir = sigm(miF-F0);
    %         mar = sigm(maF-F0);
    %         plot(hb(1),[mir,mir],[0,max(eh(:))],'color',col(ij),'linestyle','--')
    %         plot(hb(1),[mar,mar],[0,max(eh(:))],'color',col(ij),'linestyle','--')
    %     end
    may(1) = max([may(1),max(eh(:))]);
    
    
    for i=1:length(gridn)
        str{(ij-1)*2+i} = [s0{ij},': n=',num2str(gridn(i))];
    end
    
    legend(hb(1),str)
    
    hb(2) = subplot(3,1,2,'parent',hf(2),'nextplot','add');
    [eh,gy] = empiricalHist(mep(:,:,ij)',s);
    plot(hb(2),gy(:,1),eh(:,1),'color',col(ij),'linestyle','-')
    plot(hb(2),gy(:,2),eh(:,2),'color',col(ij),'linestyle','.')
    %     for i=1:length(gridn)
    %         t(i) = gy(min(find(cumsum(eh(:,i))>=1-alpha)),i);
    %         plot(hb(2),[t(i),t(i)],[0,max(eh(:))],'color',col(ij))
    %         mif = 1./K;
    %         plot(hb(2),[mif,mif],[0,max(eh(:))],'k--')
    %     end
    may(2) = max([may(2),max(eh(:))]);
    
    
    hb(3) = subplot(3,1,3,'parent',hf(2),'nextplot','add');
    [eh,gy] = empiricalHist(mEf(:,:,ij)',s);
    plot(hb(3),gy(:,1),eh(:,1),'color',col(ij),'linestyle','-')
    plot(hb(3),gy(:,2),eh(:,2),'color',col(ij),'linestyle','.')
    %     for i=1:length(gridn)
    %         t(i) = gy(min(find(cumsum(eh(:,i))>=1-alpha)),i);
    %         plot(hb(3),[t(i),t(i)],[0,max(eh(:))],'color',col(ij))
    %         mif = (ceil(gridn(i)/K)+a0)./(gridn(i)+K*a0);
    %         maf = (gridn(i)+a0)./(gridn(i)+K*a0);
    %         plot(hb(3),[mif,mif],[0,max(eh(:))],'color',col(ij),'linestyle','--')
    %         plot(hb(3),[maf,maf],[0,max(eh(:))],'color',col(ij),'linestyle','--')
    %     end
    may(3) = max([may(3),max(eh(:))]);
    
    
end

set(hb(1),'xlim',[0,1],'ylim',[0,may(1)])
title(hb(1),'BF''s histogram')
set(hb(2),'xlim',[0,1],'ylim',[0,may(2)])
title(hb(2),'max(ep)''s histogram')
set(hb(3),'xlim',[0,1],'ylim',[0,may(3)])
title('max(Ef)''s histogram')

getSubplots




