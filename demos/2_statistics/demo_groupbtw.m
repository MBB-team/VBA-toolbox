% This demo evaluates group-level BMS when asking whether different
% conditions induce different models. More precisely, the code uses MCMC
% simulations to obtain distribution of Eps, BORs and protected EPs, given
% fixed frequency profiles of tuple families f= and f~=.
% See function VBA_groupBMCbtw.m

clear all
close all

Ndummy = 128;
Nmcmc = 16;
N = 32;
r0 = [0.5;0.5]; % model frequencies under the null 
dr = [-0.5:0.05:0.5]; % change in model frequencies
ep = zeros(Nmcmc,length(dr)); % exceedance prob (f=)
pep = zeros(Nmcmc,length(dr)); % protected exceedance prob (f=)
bor = zeros(Nmcmc,length(dr)); % bayesian omnibus risk

% create dummy log-evidences
L0(:,:,1) = [ones(1,Ndummy);-ones(1,Ndummy)];
L0(:,:,2) = -L0(:,:,1);
L0 = L0 + 1e-1*randn(size(L0));

for j=1:length(dr)
    j./length(dr)
    for ii=1:Nmcmc
        % sample 2-tuple families (homogeneous vs heterogeneous conditions)
        f = VBA_random ('Categorical', r0 + [dr(j); -dr(j)], N, 1);
        
        % sample model evidences for each model, subject and condition
        L = zeros(2,N,2); % 2xNx2 (2 models; N subjects; 2 conditions).
        for i=1:N
            i1 = floor(1+rand*(Ndummy-1));
            i2 = floor(1+rand*(Ndummy-1));
            % pick a model at random
            m = VBA_random ('Categorical', [0.5; 0.5]);
            
            % assign LEVs, given the subject's class
            if f(i)==1 % homogeneous condition
                L(:,i,1) = L0(:,i1,m);
                L(:,i,2) = L0(:,i2,m);
            else % heterogeneous condition
                L(:,i,1) = L0(:,i1,m);
                L(:,i,2) = L0(:,i2,3-m);
            end
        end
        [ep(ii,j),out] = VBA_groupBMC_btwConds(L,struct('verbose',0,'DisplayWin',0));
        pep(ii,j) = out.pep;
        bor(ii,j) = out.VBA.btw.out.bor;        
    end
end


% form frequency-dependent histograms
nx = [0:0.1:1];
[ny,nx] = hist(ep,nx');
pX(:,:,1) = ny'./Nmcmc;
[ny,nx] = hist(bor,nx');
pX(:,:,2) = ny'./Nmcmc;
[ny,nx] = hist(pep,nx');
pX(:,:,3) = ny'./Nmcmc;
hf = figure('color',[1 1 1],'name','group-BMC summary statistics: MCMC histograms');
th = {'EP(1)','BOR','protected EP(1)'};
for i=1:3
    ha(i) = subplot(2,2,i,'parent',hf);
    plotGraph3D(pX(:,:,i),nx',[],ha(i))
    set(ha(i),...
        'zlim',[0,1],...
        'xtick',[0,length(dr)/2-1,length(dr)-1],...
        'xticklabel',[r0(1)+dr(1),r0(1),r0(1)+dr(end)],...
        'view',[-20,60])
    title(ha(i),th{i})
    xlabel(ha(i),'dr')
    ylabel(ha(i),th{i})
    zlabel(ha(i),'MCMC histogram')
end


% scramble frequency-dependent histogams and perform ROC analysis
clear out
ip = find(dr>0);
in = find(dr<0);
pr = 1;
y0 = VBA_vec(ep(:,in));
y1 = VBA_vec(ep(:,ip));
[h0(:,1),g0(:,1)] = VBA_empiricalDensity(y0,pr);
[h1(:,1),g1(:,1)] = VBA_empiricalDensity(y1,pr);
[proc(1),out(1),hf] = doROC(y1,y0);
set(hf,'name','ROC analysis: EP')
y0 = VBA_vec(1-bor(:,in));
y1 = VBA_vec(1-bor(:,ip));
[h0(:,2),g0(:,2)] = VBA_empiricalDensity(y0,pr);
[h1(:,2),g1(:,2)] = VBA_empiricalDensity(y1,pr);
[proc(2),out(2),hf] = doROC(y1,y0);
set(hf,'name','ROC analysis: BOR')
y0 = VBA_vec(pep(:,in));
y1 = VBA_vec(pep(:,ip));
[h0(:,3),g0(:,3)] = VBA_empiricalDensity(y0,pr);
[h1(:,3),g1(:,3)] = VBA_empiricalDensity(y1,pr);
[proc(3),out(3),hf] = doROC(y1,y0);
set(hf,'name','ROC analysis: protected EP')
hf = figure('color',[1 1 1],'name','group-BMC summary statistics: ROC analysis');
th = {'max EP','1-BOR','max protected EP'};
for i=1:3
    ha(i) = subplot(2,2,i,'parent',hf,'nextplot','add','xlim',[0,1],'ylim',[0,0.02]);
    plot(ha(i),g0(:,i),h0(:,i),'r')
    plot(ha(i),g1(:,i),h1(:,i),'g')
    xlabel(ha(i),th{i})
    ylabel(ha(i),'MCMC histogram')
    legend(ha(i),{'H0','H1'})
end


% assess self-consistency
nx = [0:0.1:1];
ip = find(dr>=0);
for ii=1:2048
    ind =  floor(1+rand(128,1)*(Nmcmc-1));
    [ny,nx] = hist(ep(ind,:),nx');
    pX2(:,:,1) = ny'./128;
    [ny,nx] = hist(pep(ind,:),nx');
    pX2(:,:,2) = ny'./128;
    for i=1:length(nx)
        pr = pX2(:,i,1);
        pr = pr./sum(pr);
        tp(i,1,ii) = sum(pr(ip));
        pr = pX2(:,i,2);
        pr = pr./sum(pr);
        tp(i,2,ii) = sum(pr(ip));
    end
end
mtp = mean(tp,3);
stp = std(tp,[],3);
hf = figure('color',[1 1 1],'name','BMC: between-condition comparison');
ha = axes('parent',hf,'nextplot','add');
plot(ha,nx,mtp(:,1),'color','r','linestyle','-','marker','.')
plot(ha,nx,mtp(:,2),'color','g','linestyle','-','marker','.')
errorbar(ha,nx,mtp(:,1),1.96*stp(:,1),'r')
errorbar(ha,nx,mtp(:,2),1.96*stp(:,2),'g')
legend(ha,{'EP','protected EP'})
plot(ha,[0,1],[0,1],'k--')
set(ha,'xlim',[0,1],'ylim',[0,1])
xlabel(ha,'sampled P(r1>r0)')
ylabel(ha,'evaluated P(r1>r0|y)')
VBA_getSubplots ();

