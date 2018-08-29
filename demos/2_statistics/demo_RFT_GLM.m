% demo for RFT simulations

clear all
close all
clc


%--- RFT: GLM analysis
L = 1e3; % size of the 1D field
n = 64; % number of GLM data samples
k = 8; % number of GLM parameters
N = 1e2; % number of Monte-Carlo simulations
v = ([1:5:50]/2.355).^2; % smoothing kernel variances
verbose = 0;
for j=1:length(v)
    j
    kernel = exp(-0.5.*([1:L]-L/2).^2/v(j));
    kernel = kernel./sum(kernel);
    for ii=1:N
        % simulate data under GLM
        e = zeros(L,n);
        for i=1:n
            e(:,i) = conv(randn(L,1),kernel,'same');
        end
        e = e./repmat(sqrt(var(e,[],2)),1,n);
        X = randn(n,k);
        b = zeros(k,L); % under the null!
        y = X*b + e';
        for i=1:n
            y(i,:) = conv(y(i,:),kernel,'same');
        end
        % apply RFT
        c = [1;zeros(k-1,1)];
        [stat,out] = RFT_GLM_contrast(X,y,c,'t',[],verbose);
        hpeak(j,ii) = length(find(out.peaks.prft<0.05))>=1;
        hclu(j,ii) = length(find(out.clusters.prft<0.05))>=1;
        f(j,ii) = out.fwhm;
        % apply Bonferroni
        p = 1-VBA_spm_Tcdf(stat,out.options.dof);
        h(j,ii) = length(find(p<0.05./L))>=1;
    end
end


hf = figure('color',[1 1 1]);
ha = subplot(2,1,1,'parent',hf,'nextplot','add');
errorbar(sqrt(v*8*log(2)),mean(hpeak'),std(hpeak')/sqrt(N),'parent',ha)
errorbar(sqrt(v*8*log(2)),mean(hclu'),std(hclu')/sqrt(N),'parent',ha,'color','r')
errorbar(sqrt(v*8*log(2)),mean(h'),std(h')/sqrt(N),'parent',ha,'color','k')
plot(ha,sqrt(v*8*log(2)),0.05*ones(size(v)),'k--')
xlabel(ha,'smoothing kernel FWHM')
ylabel(ha,'controlled FWER (test size=5%)')
legend(ha,{'RFT: peak','RFT: cluster','Bonferroni'})
ha = subplot(2,1,2,'parent',hf,'nextplot','add');
errorbar(sqrt(v*8*log(2)),mean(f'),std(f')/sqrt(N),'parent',ha)
plot(ha,sqrt(v*8*log(2)),sqrt(v*8*log(2)),'k--')
xlabel(ha,'smoothing kernel FWHM')
ylabel(ha,'RFT-FWHM')
VBA_getSubplots ();