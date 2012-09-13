function [out] = group_level_analysis(lme,type,model_names,partition,alpha0, Nsamp)
% Performs Bayesian model comparison at the group level
% Model comparison can be performed under the 2 following assumptions
% - FFX : All subjects use the same model
% - RFX : There is a probability at the population level on each model
% being used

% INPUT
% - lme : [Nsubjects x Nmodels] matrix of log evidences
% - type : 'RFX' or 'FFX', type of inference
% - model_names : cell of string of model names (can be left empty)
% OUTPUT
% - out : structure
%   if type is 'RFX':
%       .exp_r   : [1 x  Nmodels] expectation of the posterior p(r|y)
%       .xp      : exceedance probabilities
%       .r_samp  : [Nsamp x Nmodels] matrix of samples from posterior
%       .g_post  : [Nsubject x Nmodels] matrix of posterior probabilities with
%                   g_post(i,k) being post prob that subj i used model k
%   if type is 'FFX':
%       .pp      : posterior probability on models
%       .lme_all : for each model, summed log evidence over subjects


disp('Starting group level analysis...')

Nm = size(lme,2); % number of models

try 
    partition;
catch
    partition = cell(1,Nm);
    for i = 1:Nm
        partition{i}=i;
    end 
end

[lme] = models2families(lme,partition);
Nm = size(lme,2); % number of models

% ------------------------------------ RANDOM EFFECTS ANALYSIS

if isequal(type,'RFX') % Random effect analysis
    
    try alpha0;
    catch; alpha0 = ones(1,Nm); % prior count on models
    end
    
    try Nsamp;
    catch; Nsamp = 1e4; % number of gibbs samples
    end
    
    [exp_r,xp,r_samp,g_post] = spm_BMS_gibbs (lme, alpha0, Nsamp);
    out.exp_r = exp_r;
    out.xp = xp;
    try  out.rsamp = rsamp;catch;end
    try  out.g_post = g_post;catch;end
    
    
    
    f= figure;
    Pos = get(f,'Position');
    Pos(3:4) =  [1200,400];
    set(f, 'Position', Pos); % [left bottom width height]
    set(f,'name','Random Effect Analysis','numbertitle','off')
    
    
    %-- plotting log-evidences
    subplot(1,3,1);
    b= bar(lme');
    set(b,'barwidth',1)
    
    
    title('Log evidences per model', ...
        'FontSize', 10, ...
        'FontWeight', 'bold')
    
    ylabel('Log evidences', ...
        'FontSize', 10, ...
        'FontWeight', 'bold')
    
    % if exist('model_names','var')
    try
        set(gca, 'XTickLabel',model_names)
        th = rotateticklabel(gca,45);
    catch
        xlabel('Models', ...
            'FontSize', 10, ...
            'FontWeight', 'bold')
    end
    
    %-- plotting exceedance probabilities
    subplot(1,3,2)
    h1=bar(xp);
    set(h1,'facecolor','black')
    title('Exceedance probabilities', ...
        'FontSize', 10, ...
        'FontWeight', 'bold')
    
    ylabel('Exceedance probability', ...
        'FontSize', 10, ...
        'FontWeight', 'bold')
    
    if exist('model_names','var')
        try
            set(gca, 'XTickLabel',model_names)
            th = rotateticklabel(gca,45);
        catch; 
        xlabel('Models', ...
            'FontSize', 10, ...
            'FontWeight', 'bold')
        end
    end
    
    
    %-- plotting expectation of the posterior p(r|y)
    subplot(1,3,3)
    h2=bar(exp_r);
    set(h2,'facecolor','black')
    title('Expected model frequencies', ...
        'FontSize', 10, ...
        'FontWeight', 'bold')
    
    ylabel('p(r_i|y)', ...
        'FontSize', 10, ...
        'FontWeight', 'bold')
    
    if exist('model_names','var')
        try
            set(gca, 'XTickLabel',model_names)
            th = rotateticklabel(gca,45);
        catch;end;
    else
        xlabel('Models', ...
            'FontSize', 10, ...
            'FontWeight', 'bold')
    end
    disp('Done!')
    
    % ------------------------------------ FIXED EFFECTS ANALYSIS
    
elseif isequal(type,'FFX') % fixed effects analysis
    
    lme_all = sum(lme,1);
    pp = exp(lme_all-max(lme_all))./sum(exp(lme_all-max(lme_all)));
    out.pp = pp;
    out.lme_all = lme_all;
    
 
    
    f= figure;
    Pos = get(f,'Position');
    Pos(3:4) =  [800,400];
    set(f, 'Position', Pos); % [left bottom width height]
    set(f,'name','Fixed Effect Analysis','numbertitle','off')
    
    
    %-- plotting log-evidences for each subject and model
    subplot(1,3,1)
    b= bar(lme');
    set(b,'barwidth',1)
    
    title('Log evidences per model', ...
        'FontSize', 10, ...
        'FontWeight', 'bold')
    
    ylabel('Log evidences', ...
        'FontSize', 10, ...
        'FontWeight', 'bold')
    
    %      l =  get(gca,'YLim');
    %      l(1) = min(min(lme));
    %      l(2) = max(max(lme));
    %      set(gca,'YLim',l)
    %  %    get(gca,'YTick')
    
    
    if exist('model_names','var')
        try
            set(gca,'Layer','top','XTickLabel',model_names)
            th = rotateticklabel(gca,45);
        catch;end
    else
        set(gca,'YLim',[min(min(lme)) max(max(lme))])
        xlabel('Models', ...
            'FontSize', 10, ...
            'FontWeight', 'bold')
    end
    
    
    
    %-- plotting log-evidences for each model
    
    subplot(1,3,2)
    h1=bar(lme_all);
    set(h1,'facecolor','black')
    title('Model Evidence', ...
        'FontSize', 10, ...
        'FontWeight', 'bold')
    
    ylabel('Log Evidence', ...
        'FontSize', 10, ...
        'FontWeight', 'bold')
    
    
    % %          l =  get(gca,'YLim');
    %      l(1) = min(min(lme_all));
    %      l(2) = max(max(lme_all));
    %      set(gca,'YLim',l)
    % %     get(gca,'YTick')
    
    %if exist('model_names','var')
    try
        set(gca, 'XTickLabel',model_names)
        th = rotateticklabel(gca,45);
    catch%else
        set(gca,'YLim',[min(min(lme_all)) max(max(lme_all))])
        xlabel('Models', ...
            'FontSize', 10, ...
            'FontWeight', 'bold')
    end
    
    disp('Done!')
    
    %-- plotting posterior on each model
    subplot(1,3,3)
    h1=bar(pp);
    set(h1,'facecolor','black')
    title('Posterior on model', ...
        'FontSize', 10, ...
        'FontWeight', 'bold')
    
    ylabel('Posterior probability', ...
        'FontSize', 10, ...
        'FontWeight', 'bold')
    
    if exist('model_names','var')
    try
        set(gca, 'XTickLabel',model_names)
        th = rotateticklabel(gca,45);
    catch;end;
    else
        xlabel('Models', ...
            'FontSize', 10, ...
            'FontWeight', 'bold')
    end
    disp('Done!')
    
    
else
    % Do nothing
    disp('Could not perform group level analysis, check your inputs!')
end


end

function th=rotateticklabel(h,rot)
%ROTATETICKLABEL rotates tick labels
%   Written Oct 14, 2005 by Andy Bliss
%   Copyright 2005 by Andy Bliss
%set the default rotation if user doesn't specify
if nargin==1
    rot=90;
end
%make sure the rotation is in the range 0:360 (brute force method)
while rot>360
    rot=rot-360;end
while rot<0
    rot=rot+360;end
a=get(h,'XTickLabel');%get current tick labels
set(h,'XTickLabel',[]);%erase current tick labels from figure
b=get(h,'XTick');%get tick label positions
c=get(h,'YTick');%make new tick labels

if rot<180
    th=text(b,repmat(c(1)-.1*(c(2)-c(1)),length(b),1),a,'HorizontalAlignment','right','rotation',rot);
else
    th=text(b,repmat(c(1)-.1*(c(2)-c(1)),length(b),1),a,'HorizontalAlignment','left','rotation',rot);
end
end