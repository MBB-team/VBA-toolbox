% This demo inspects the properties of frequentist statistical tests (e.g.
% power) as a function of sample size.
% In brief, we test for the mean of a sample with varying sample size,
% under the null (mean=0) and under the alternative hypotheses (mean=1,
% here). We then compute the FPR, FNR, TPR and TNR, as well as
% sensitivity/specificity of the test.
% NB: we also vary the hidden effect size from a very small
% ('uninteresting') to a 'normal' magnitude. We then look at summary
% statistics that correct for 'uninteresting' findings...



gridn = 2.^[1:8];
gridr = [1e-2;1e-1;1;1e1];
N = 2^11;

FP = zeros(length(gridn),length(gridr));
FN = FP;
TP = FP;
TN = FP;

alpha = 0.05;


for i=1:N
    
    for j=1:length(gridn)
        
        n =gridn(j);
        X = ones(n,1);
        e = randn(n,1);
%         e = e - mean(e);
        
        for k=1:length(gridr)
            
            r = gridr(k);
            y1 = r*X + e;
            y0 = e;
            
            [p1] = GLM_contrast(X,y1,1,'t');
            [p0] = GLM_contrast(X,y0,1,'t');
            
            if p1<=alpha
                TP(j,k) = TP(j,k) +1;
            else
                FN(j,k) = FN(j,k) +1;
            end
            
            if p0<=alpha
                FP(j,k) = FP(j,k) +1;
            else
                TN(j,k) = TN(j,k) +1;
            end
            
        end
        
    end
    
end
sen = TP./(TP+FN);
spe = TN./(TN+FP);
PPV = TP./(TP+FP);
NPV = TN./(TN+FN);

FPR = FP/N;
TPR = TP/N;
FNR = FN/N;
TNR = TN/N;

hf = figure('color',[1 1 1]);
for k=1:length(gridr)
    
    ha = subplot(...
        length(gridr),2,2*(k-1)+1,...
        'parent',hf,...
        'nextplot','add',...
        'ylim',[0 1]);
    plot(ha,FPR(:,k),'b')
    plot(ha,FNR(:,k),'g');
    plot(ha,TPR(:,k),'r')
    plot(ha,TNR(:,k),'m')
    legend({'FPR = Type I error rate (controlled)','FNR = Type II error rate (1-power)','TPR','TNR'})
    set(ha,'xtick',1:length(gridn),'xticklabel',gridn)
    xlabel(ha,'sample size')
    title(ha,['effect size:',num2str(gridr(k))])
    
    ha = subplot(...
        length(gridr),2,2*(k-1)+2,...
        'parent',hf,...
        'nextplot','add',...
        'ylim',[0 1]);
    plot(ha,sen(:,k),'b')
    plot(ha,spe(:,k),'g');
    plot(ha,PPV(:,k),'r');
    plot(ha,NPV(:,k),'m');
    legend({'sensitivity','specificity','PPV','NPV'})
    set(ha,'xtick',1:length(gridn),'xticklabel',gridn)
    xlabel(ha,'sample size')
    title(ha,['effect size:',num2str(gridr(k))])
    
end

% Now reconsider true effect sizes below a given a threhold (rt) as being
% 'uninteresting'. This means that for r<rt, TP are actually FP, and FN are
% actually TN. One can then average across effect sizes, and produce
% 'reasonable' summary statistics, as follows:

hf2 = figure(...
        'color',[1 1 1],...
        'name','The fallacy of classical statistics!');

for k=1:length(gridr)
    
    rt = gridr(k);
    
    inr = find(gridr>=rt); % 'interesting' effects
    outr = find(gridr<rt); % 'un-interesting' effects
    
    TP2 = sum(TP(:,inr),2);
    FN2 = sum(FN(:,inr),2);
    FP2 = sum(FP,2) + sum(TP(:,outr),2);
    TN2 = sum(TN,2) + sum(FN(:,outr),2);
    
    Nout = TN2+FP2;% number of corrected true null simulations
    Nin = TP2+FN2; % number of corrected true alternative simulations
    
    sen2 = TP2./Nin;
    spe2 = TN2./Nout;
    PPV2 = TP2./(TP2+FP2);
    NPV2 = TN2./(TN2+FN2);
    
    
    FPR2 = FP2./Nout;
    TPR2 = TP2./Nin;
    FNR2 = FN2./Nin;
    TNR2 = TN2./Nout;
    
    
    ha = subplot(...
        length(gridr),2,2*(k-1)+1,...
        'parent',hf2,...
        'nextplot','add',...
        'ylim',[0 1]);
    plot(ha,FPR2,'b')
    plot(ha,FNR2,'g');
    plot(ha,TPR2,'r')
    plot(ha,TNR2,'m')
    legend({'FPR = Type I error rate (controlled)','FNR = Type II error rate (1-power)','TPR','TNR'})
    set(ha,'xtick',1:length(gridn),'xticklabel',gridn)
    xlabel(ha,'sample size')
    title(ha,['effect size threshold:',num2str(rt)])
    
    ha = subplot(...
        length(gridr),2,2*(k-1)+2,...
        'parent',hf2,...
        'nextplot','add',...
        'ylim',[0 1]);
    plot(ha,sen2,'b')
    plot(ha,spe2,'g');
    plot(ha,PPV2,'r');
    plot(ha,NPV2,'m');
    legend({'sensitivity','specificity','PPV','NPV'})
    set(ha,'xtick',1:length(gridn),'xticklabel',gridn)
    xlabel(ha,'sample size')
    title(ha,['effect size threshold:',num2str(rt)])
    
end

