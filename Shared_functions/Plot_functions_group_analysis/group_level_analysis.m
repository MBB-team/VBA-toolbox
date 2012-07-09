% Aim of the function
% - perform RFX and FFX analysis
% - Graphic display of the results




function [] = group_level_analysis(lme)
%------------------------------------------------
% Performing group level analysis
%------------------------------------------------

disp('Starting group level analysis...')

Ns = size(lme,1); % number of subjects
Nm = size(lme,2); % number of models

% Random effect analysismatlab

alpha0 = ones(1,Nm); % prior count on models
Nsamp = 1e4;
[exp_r,xp,r_samp,g_post] = spm_BMS_gibbs (lme, alpha0, Nsamp);

%%
%------------------------------------------------
% Display results
%------------------------------------------------
% exp_r   - [1 x  Nk] expectation of the posterior p(r|y)
% xp      - exceedance probabilities
% r_samp  - [Nsamp x Nk] matrix of samples from posterior
% g_post  - [Ni x Nk] matrix of posterior probabilities with 
%           g_post(i,k) being post prob that subj i used model k

f= figure;
Pos = get(f,'Position');
Pos(3:4) =  [1200,400];
set(f, 'Position', Pos); % [left bottom width height]
set(f,'name','Random Effect Analysis','numbertitle','off')


%-- plotting log-evidences
s = subplot(1,3,1)
b= bar(lme');
set(b,'barwidth',1)


title('Log evidences per model', ...
'FontSize', 10, ...
'FontWeight', 'bold')
xlabel( 'Models', ...
'FontSize', 10, ...
'FontWeight', 'bold')
ylabel('Log evidences', ...
'FontSize', 10, ...
'FontWeight', 'bold')


%-- plotting exceedance probabilities
subplot(1,3,2)
h1=bar(xp);
set(h1,'facecolor','black')
title('Exceedance probabilities', ...
'FontSize', 10, ...
'FontWeight', 'bold')
xlabel('Models', ...
'FontSize', 10, ...
'FontWeight', 'bold')
ylabel('Exceedance probability', ...
'FontSize', 10, ...
'FontWeight', 'bold')

%-- plotting expectation of the posterior p(r|y)
subplot(1,3,3)
h2=bar(exp_r);
set(h2,'facecolor','black')
title('Expectation of the posterior p(r|y)', ...
'FontSize', 10, ...
'FontWeight', 'bold')
xlabel('Models', ...
'FontSize', 10, ...
'FontWeight', 'bold')
ylabel('p(r_i|y)', ...
'FontSize', 10, ...
'FontWeight', 'bold')

%%



f= figure;
Pos = get(f,'Position');
Pos(3:4) =  [800,400];
set(f, 'Position', Pos); % [left bottom width height]
set(f,'name','Fixed Effect Analysis','numbertitle','off')


%-- plotting log-evidences
s = subplot(1,2,1)
b= bar(lme');
set(b,'barwidth',1)

title('Log evidences per model', ...
'FontSize', 10, ...
'FontWeight', 'bold')
xlabel( 'Models', ...
'FontSize', 10, ...
'FontWeight', 'bold')
ylabel('Log evidences', ...
'FontSize', 10, ...
'FontWeight', 'bold')


%-- Log Evidences
subplot(1,2,2)
h1=bar(sum(lme,1));
set(h1,'facecolor','black')
title('Log Evidence', ...
'FontSize', 10, ...
'FontWeight', 'bold')
xlabel('Models', ...
'FontSize', 10, ...
'FontWeight', 'bold')
ylabel('Log Evidence', ...
'FontSize', 10, ...
'FontWeight', 'bold')


disp('Done!')

end

