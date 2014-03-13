% Polya's urn

% close all
clear all
clc

flag = 'cep'; % 'Ef'

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
    if isequal(flag,'Ef')
        ylabel(ha(i,2),'estimated freq')
    elseif isequal(flag,'cep')
        ylabel(ha(i,2),'corrected EP')
    end
    mif = 1./K;
    maf = (gridn(i)+a0)./(gridn(i)+K*a0);
    plot(ha(i,2),[0,1],[mif,mif],'b--')
    plot(ha(i,2),[0,1],[maf,maf],'b--')
    plot(ha(i,2),[mir,mir],[0,1],'b--')
    plot(ha(i,2),[mar,mar],[0,1],'b--')
    
    ha(i,3) = subplot(length(gridn),3,(i-1)*3+3,'parent',hf,'nextplot','add','xlim',[0,1],'ylim',[0,1]);
    title(ha(i,3),str{i})
    if isequal(flag,'Ef')
        xlabel(ha(i,3),'estimated freq')
    elseif isequal(flag,'cep')
        xlabel(ha(i,3),'corrected EP')
    end
    ylabel(ha(i,3),'exceedance prob')
    plot(ha(i,3),[0,1],[1,1]./K,'b--')
    plot(ha(i,3),[mif,mif],[0,1],'b--')
    plot(ha(i,3),[maf,maf],[0,1],'b--')
end

FPR = zeros(length(gridn),2);
r = zeros(length(gridn),N,2);
mep = zeros(length(gridn),N,2);
cep = zeros(length(gridn),N,2);
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
        cep(i,j,1) = mep(i,j,1).*r(i,j,1);
        mEf(i,j,1) = max(Ef);
        % invert model under the alternative
        f2 = VBA_sample('dirichlet',struct('d',a0*ones(K,1)),1,0);
        m = sampleFromArbitraryP(f2,eye(length(f)),n)';
        [F,F0,a] = PolyaLEV(m,a0);
        ep = VBA_ExceedanceProb(a,[],'dirichlet',0);
        Ef = a./sum(a);
        r(i,j,2) = sigm(F-F0);
        mep(i,j,2) = max(ep);
        cep(i,j,2) = mep(i,j,2).*r(i,j,2);
        mEf(i,j,2) = max(Ef);
        % display statistics
        plot(ha(i,1),r(i,j,1),mep(i,j,1),'r.')
        plot(ha(i,1),r(i,j,2),mep(i,j,2),'g.')
        if isequal(flag,'Ef')
            plot(ha(i,2),r(i,j,1),mEf(i,j,1),'r.')
            plot(ha(i,3),mEf(i,j,1),mep(i,j,1),'r.')
            plot(ha(i,2),r(i,j,2),mEf(i,j,2),'g.')
            plot(ha(i,3),mEf(i,j,2),mep(i,j,2),'g.')
        elseif isequal(flag,'cep')
            plot(ha(i,2),r(i,j,1),cep(i,j,1),'r.')
            plot(ha(i,3),cep(i,j,1),mep(i,j,1),'r.')
            plot(ha(i,2),r(i,j,2),cep(i,j,2),'g.')
            plot(ha(i,3),cep(i,j,2),mep(i,j,2),'g.')
        end
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
    may(1) = max([may(1),max(eh(:))]);
    
    
    for i=1:length(gridn)
        str{(ij-1)*2+i} = [s0{ij},': n=',num2str(gridn(i))];
    end
    legend(hb(1),str)
    hb(2) = subplot(3,1,2,'parent',hf(2),'nextplot','add');
    [eh,gy] = empiricalHist(mep(:,:,ij)',s);
    plot(hb(2),gy(:,1),eh(:,1),'color',col(ij),'linestyle','-')
    plot(hb(2),gy(:,2),eh(:,2),'color',col(ij),'linestyle','.')
    may(2) = max([may(2),max(eh(:))]);
    
    
    
    if isequal(flag,'Ef')
        [eh,gy] = empiricalHist(mEf(:,:,ij)',s);
    elseif isequal(flag,'cep')
        [eh,gy] = empiricalHist(cep(:,:,ij)',s);
    end
    hb(3) = subplot(3,1,3,'parent',hf(2),'nextplot','add');
    plot(hb(3),gy(:,1),eh(:,1),'color',col(ij),'linestyle','-')
    plot(hb(3),gy(:,2),eh(:,2),'color',col(ij),'linestyle','.')
    may(3) = max([may(3),max(eh(:))]);
    
    
end
set(hb(1),'xlim',[0,1],'ylim',[0,may(1)])
title(hb(1),'BF''s histogram')
set(hb(2),'xlim',[0,1],'ylim',[0,may(2)])
title(hb(2),'max(ep)''s histogram')
set(hb(3),'xlim',[0,1],'ylim',[0,may(3)])
if isequal(flag,'Ef')
    title(hb(3),'max(Ef)''s histogram')
elseif isequal(flag,'cep')
    title(hb(3),'max(cEf)''s histogram')
end

for i=1:length(gridn)
    
    xp = mep(i,:,2);
    xn = mep(i,:,1);
    [p,out] = doROC(xp,xn);
    set(gcf,'name',['max(ep) ROC analysis: n=',num2str(gridn(i))])
    
    if isequal(flag,'Ef')
        xp = mEf(i,:,2);
        xn = mEf(i,:,1);
    elseif isequal(flag,'cep')
        xp = cep(i,:,2);
        xn = cep(i,:,1);
    end
    [p,out] = doROC(xp,xn);
    if isequal(flag,'Ef')
        set(gcf,'name',['max(Ef) ROC analysis: n=',num2str(gridn(i))])
    elseif isequal(flag,'cep')
        set(gcf,'name',['max(cEP) ROC analysis: n=',num2str(gridn(i))])
    end
    xp = r(i,:,2);
    xn = r(i,:,1);
    [p,out] = doROC(xp,xn);
    set(gcf,'name',['BF ROC analysis: n=',num2str(gridn(i))])
end

getSubplots

