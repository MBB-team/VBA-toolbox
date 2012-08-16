function [] = VBA_Redisplay_modulation(LogEv,mod)
% This function provides a graphical presentation of the model comparison
% of nested models in which modulatory parameters are
% activated/deactivated.
% INPUT
% - LogEv : matrix of log evidences for each model and subject (Nmodels*Nsubjects)
% - inv_order : matrix specifying the status (activated/deactivated) of
% each parametric modulator for each models
% - modulator names : names of each of the parametric modulators ordered as
% in inv_order
% OUTPUT : none

modulator_names = mod.modulator_names;
inv_order = mod.inv_order;
inv_order_gm = mod.inv_order_gm;
N_gm = size(inv_order_gm,2);
Nmodels = size(inv_order_gm,1);


Nsubject = size(LogEv,2); % number of subjects considered

%------------------- Case single subject ---------------------


if Nsubject == 1 % Case single subject


Nplots = 4;

% if isempty(modulator_names)
%     modulator_names = cell(1,N_p);
% end

pp = exp(LogEv-max(LogEv));
pp = pp./sum(pp);

figure

%----------
subplot(1,Nplots,1)
imagesc(-inv_order_gm)
colormap(gray)

if ~isempty(modulator_names)

N_m = length(modulator_names); % number of modulators
set(gca, 'XTickLabel',modulator_names) % setting label names
set(gca,'XLim',[0.5 N_m+0.5]) 
set(gca,'XTick',[1:N_m]) % setting label positions

set(gca,'XTickLabel',[]);%erase current tick labels from figure
a=get(gca,'XTickLabel');%make new tick labels

c=get(gca,'YTick')+N_gm-1;%make new tick labels
b=get(gca,'XTick');%get tick label positions
text(b,repmat(N_gm+1,N_m,1),modulator_names,'HorizontalAlignment','right','rotation',45); % rotating labels 
 
else
    
    xlabel('Modulation parameters', ...
        'FontSize', 10, ...
        'FontWeight', 'bold')
end


ylabel('Models')

%----------
subplot(1,Nplots,2)
barh((LogEv-max(LogEv)),0.5)
axis([min((LogEv-max(LogEv))) max((LogEv-max(LogEv))) 1-0.5 Nmodels+0.5])
set(gca, 'YDir', 'reverse')
xlabel('Log Evidence')

%----------
subplot(1,Nplots,3)
barh(pp,0.5)
axis([0 max(pp) 1-0.5 Nmodels+0.5])
set(gca, 'YDir', 'reverse')
xlabel('Posterior probability')

% marginalized
ppm = zeros(1,N_gm);
for i_g = 1:N_gm
  ppm(i_g) = sum(pp(find(inv_order(:,i_g))));
end

%----------
subplot(1,Nplots,4)
bar(ppm,0.5)
axis([ 1-0.5 N_gm+0.5 0 max(ppm)])
xlabel('Marginal probabilities')

%------------------- Case multiple subjects ---------------------

else % Case multiple subject
    
Nplots = 4;


% performing random effects analysis
[exp_r,xp,r_samp,g_post] = spm_BMS_gibbs (LogEv, ones(1,Nmodels))



figure

%----------
subplot(1,Nplots,1)

imagesc(-inv_order_gm)
colormap(gray)

if ~isempty(modulator_names)

N_m = length(modulator_names); % number of modulators
set(gca, 'XTickLabel',modulator_names) % setting label names
set(gca,'XLim',[0.5 N_m+0.5]) 
set(gca,'XTick',[1:N_m]) % setting label positions

set(gca,'XTickLabel',[]);%erase current tick labels from figure
a=get(gca,'XTickLabel');%make new tick labels

c=get(gca,'YTick')+N_gm-1;%make new tick labels
b=get(gca,'XTick');%get tick label positions
text(b,repmat(N_gm+1,N_m,1),modulator_names,'HorizontalAlignment','right','rotation',45); % rotating labels 
 
else
    
    xlabel('Modulation parameters', ...
        'FontSize', 10, ...
        'FontWeight', 'bold')
end


ylabel('Models')

%----------
subplot(1,Nplots,2)
barh((LogEv-max(max(LogEv))),0.5)
axis([min((LogEv-max(LogEv))) max((LogEv-max(LogEv))) 1-0.5 Nmodels+0.5])
set(gca, 'YDir', 'reverse')
xlabel('Log Evidence')

%----------
subplot(1,Nplots,3)
barh(xp,0.5)
axis([0 max(xp) 1-0.5 Nmodels+0.5])
set(gca, 'YDir', 'reverse')
xlabel('Posterior probability')

% marginalized
xpm = zeros(1,N_gm);
for i_g = 1:N_gm
  xpm(i_g) = sum(xp(find(inv_order(:,i_g))));
end

%----------
subplot(1,Nplots,4)
bar(xpm,0.5)
axis([ 1-0.5 N_gm+0.5 0 max(ppm)])
xlabel('Marginal probabilities')    
    
    
    
end

end

