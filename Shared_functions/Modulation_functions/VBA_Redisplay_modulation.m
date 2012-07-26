function [] = VBA_Redisplay_modulation(LogEv,inv_order,modulator_names)

Nsubject = size(LogEv,2); % number of subjects considered

if Nsubject == 1 % Case single subject
N_p = size(inv_order,1);


Nmodels = N_p;

if isempty(modulator_names)
    modulator_names = cell(1,N_p);
end

pp = exp(LogEv-max(LogEv));
pp = pp./sum(pp);

figure

subplot(1,3,1)
imagesc(-inv_order)
colormap(gray)

if ~isempty(modulator_names)

N_m = size(modulator_names,1); % number of modulators
set(gca, 'XTickLabel',modulator_names) % setting label names
set(gca,'XLim',[0.5 N_m+0.5]) 
set(gca,'XTick',[1:N_m]) % setting label positions

set(gca,'XTickLabel',[]);%erase current tick labels from figure
a=get(gca,'XTickLabel');%make new tick labels

c=get(gca,'YTick')+N_p-1;%make new tick labels
b=get(gca,'XTick');%get tick label positions
text(b,repmat(N_p+1,N_m,1),modulator_names,'HorizontalAlignment','right','rotation',45); % rotating labels 
 
else
    xlabel('Modulation parameters', ...
        'FontSize', 10, ...
        'FontWeight', 'bold')
end


ylabel('Models')

subplot(1,3,2)
barh((LogEv-max(LogEv)),0.5)
axis([min((LogEv-max(LogEv))) max((LogEv-max(LogEv))) 1-0.5 Nmodels+0.5])
set(gca, 'YDir', 'reverse')
xlabel('Log Evidence')

subplot(1,3,3)
barh(pp,0.5)
axis([0 max(pp) 1-0.5 Nmodels+0.5])
set(gca, 'YDir', 'reverse')
xlabel('Posterior probability')


else % Case multiple subject
    
    N_p = size(inv_order,1);
Nmodels = 2^N_p;

pp = exp(LogEv-max(LogEv));
pp = pp./sum(pp);

figure
subplot(1,3,1)
imagesc(-inv_order)
colormap(gray)
xlabel('Modulation parameters')
ylabel('Models')


S_LogEv = sum(LogEv,2); % sum of log-evidence
subplot(1,3,2)
barh((S_LogEv-max(S_LogEv)),0.5)
axis([min((S_LogEv-max(S_LogEv))) max((S_LogEv-max(S_LogEv))) 1-0.5 Nmodels+0.5])
set(gca, 'YDir', 'reverse')
xlabel('Log Evidence')

% performing fixed effects analysis
[exp_r,xp,r_samp,g_post] = spm_BMS_gibbs (LogEv, ones(1,Nmodels))

subplot(1,3,3)
barh(xp,0.5)
axis([0 1 1-0.5 Nmodels+0.5])
set(gca, 'YDir', 'reverse')
xlabel('Exceedance probability')

    
    
end

end

