% Polya's urn

% close all
clear all
clc

gridn = 2.^[1:2:10];
N = 256;

a0 = 1;

hf = figure('color',[1 1 1 ]);
for i=1:length(gridn)
    ha(i) = subplot(length(gridn),1,i,'parent',hf,'nextplot','add','xlim',[0,1],'ylim',[0,1]);
    title(ha(i),['n=',num2str(gridn(i))])
    plot(ha(i),[0,1],[0.05,0.05])
    plot(ha(i),[0,1],[0.95,0.95])
end

FPR = zeros(length(gridn),2);
for i=1:length(gridn)
    n = gridn(i);
    for j=1:N
        m = randn(1,n)>0;
        m = [m;1-m];
        [F,F0,a] = PolyaLEV(m,a0);
        ep = VBA_ExceedanceProb(a,[],'dirichlet',0);
        r = sigm(F-F0);
        plot(ha(i),r,ep,'r.')
        drawnow
        if max(ep)>=0.95
            FPR(i,1) = FPR(i,1) +1;
        end
        if r>=0.95
            FPR(i,2) = FPR(i,2) +1;
        end
    end
    
end



getSubplots