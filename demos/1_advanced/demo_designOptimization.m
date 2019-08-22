function demo_designOptimization ()
% // VBA toolbox //////////////////////////////////////////////////////////
%
% [posterior, out] = demo_designOptimization ()
% demo of off-line and online design optimisation
%
% This demo simulates a psychophysics paradigm similar to a signal
% detection task, whereby the detection probability is a sigmoidal function
% of the stimulus contrast (which is the design control variable).
% Our goal is to estimate  the inflexion point (detection threshold) and the
% sigmoid steepness (d prime) or the response function.
% In order to provide the most efficient estimate of these model parameters,
% we will show how to use offline (before the experiment) and online 
% (during the experiment) design optimization to decide given trial-by-trial 
% subjects' binary choice data (seen/unseen) which stimulus to show next.
%
% /////////////////////////////////////////////////////////////////////////

%% Global values
% =========================================================================

% number of trials for the respective simulations
N = 50; 

% true parameter values for the simulations 
phi = [- .5; 2.5]; % simulated parameters: [inflexion point, log-slope]

% range of potential simuli
uRange = linspace(- 1, 1, N);

% prepare display
VBA_figure('Name', 'demo_designOptimisation');

%% Define the model
% =========================================================================

% observation function 
% -------------------------------------------------------------------------
function [gx, dgdx, dgdp] = g_psychometric (~, phi, u, ~)
    % tip: VBA_sigmoid returns derivatives wrt parameters in alphabetical
    % order, ie. center -> slope.
    [gx, ~, dgdp] = VBA_sigmoid (u,...
    'center', phi(1), ...
    'slope', exp (phi(2)) ...
    );
    dgdp(2,:) = dgdp(2,:) * exp (phi(2));
    dgdx = [];
end

% dimensions
% -------------------------------------------------------------------------
dim.n_phi = 2;
dim.n_t = 1;
dim.p = N;

% general options
% -------------------------------------------------------------------------
% binary observations  
options.sources.type = 1; 

% no display
options.DisplayWin = 0;
options.verbose = 0;

%% 1) No optimisation
% =========================================================================
% here, we'll use the most naive approach, a full swipe of all possible 
% stimulus intensities

% experimental design
% -------------------------------------------------------------------------
u = uRange';

% simulate responses
% -------------------------------------------------------------------------
y = VBA_simulate (1,[],@g_psychometric,[],phi,u,[],[],options); 

% estimate parameters
% -------------------------------------------------------------------------
posterior_naive = VBA_NLStateSpaceModel (y,u,[],@g_psychometric,dim,options);

% display results
% -------------------------------------------------------------------------
plot_design(1, 'no optimisation', u, y, posterior_naive);

%% 2) Offline optimisation
% =========================================================================
% Here, we will try to find a better design before running the experiment

% experimental design
% -------------------------------------------------------------------------
% number of designs to try
nAttempts = 1e4;

% initialization
fprintf('Offline optimisation: optimizing (  0%%)');
efficiency_offOpt = - Inf;
efficiencyDesign = nan(1, nAttempts);
keepDesign = [1];

% loop over designs
for attempt = 1 : nAttempts
    
    % draw random design
    u_attempt = uRange(randi (numel (uRange), 1, N))';
        
    % estimate efficiency
    efficiencyDesign(attempt) = VBA_designEfficiency([],@g_psychometric,dim,options,u_attempt,'parameters');

    % if better, store and display
    if efficiencyDesign(attempt) > efficiency_offOpt 
        efficiency_offOpt = efficiencyDesign(attempt);
        u = sort(u_attempt);
        keepDesign(end+1) = attempt;
    end
    
    if efficiencyDesign(attempt) > efficiency_offOpt || mod(attempt, 50) == 0
        plot_design(2, 'offline optimisation', u, y, [], [], efficiencyDesign(keepDesign));
    end
    
    % progress bar
    fprintf('\b\b\b\b\b%3d%%)', round(100* attempt / nAttempts));
end

% simulate responses
% -------------------------------------------------------------------------
y = VBA_simulate (1,[],@g_psychometric,[],phi,u,[],[],options); 

% estimate parameters
% -------------------------------------------------------------------------
posterior_offline = VBA_NLStateSpaceModel (y,u,[],@g_psychometric,dim,options);

% display results
% -------------------------------------------------------------------------
plot_design(2, 'offline optimisation', u, y, posterior_offline);

%% 3) Online optimisation
% =========================================================================
% Here, we will optimize the design during the experiment by taking into
% account trial-by-trial subject's responses to adaptively select the next
% best stimulus to present

% initialization
opt = options;
u = nan(N, 1);
efficiencyInput = nan(1, length (uRange));
efficiencyDesign = nan(1, length (uRange));

% run experiment
for t = 1 : N
    
    % extend design
    % ---------------------------------------------------------------------
    % start from current posterior belief
    try
        opt.priors = posterior_online;
    end
        
    % compute efficiency of potential stimuli
    dim.p = 1;
    for i = 1 : length (uRange)
        efficiencyInput(i) = VBA_designEfficiency([],@g_psychometric,dim,opt,uRange(i),'parameters');
    end
    
    % find best next stimulus to present
    [efficiencyDesign(t), idxMaxEff] = max (efficiencyInput);
    u(t) = uRange(idxMaxEff);

    % simulate 1 responses
    % ---------------------------------------------------------------------
    y(t) = VBA_simulate (1,[],@g_psychometric,[],phi,u(t),Inf,[],options);
    
    % estimate parameters given data acquired so far
    % ---------------------------------------------------------------------
    dim.p = t;
    posterior_online = VBA_NLStateSpaceModel(y(1:t),u(1:t),[],@g_psychometric,dim,options);

    % display
    % ---------------------------------------------------------------------
    plot_design(3, 'online optimisation', u, y, posterior_online, efficiencyInput, efficiencyDesign);

end

% display
% ---------------------------------------------------------------------  
plot_design(3, 'online optimisation', u, y, posterior_online);

%% show results
% =========================================================================

fprintf('\nSimulation results:\n');

disp (table ( ...
    phi, ...
    posterior_naive.muPhi, ...
    posterior_offline.muPhi, ...
    posterior_online.muPhi, ...
    'RowNames', {'center','slope'}, ...
    'VariableNames',{'true','naive','offline','online'}));

%% ########################################################################
%  display subfunction
%  ########################################################################

function plot_design(idx, titleTxt, u, y, posterior, uEfficiency, dEfficiency)
        
    % jitter for data display
    persistent jitter;
    if isempty(jitter)
        jitter = 0.1 * (rand(N,1)-0.5);
    end

    % experimental design 
    subplot(4,3,idx)
    
    if nargin > 5 && ~ isempty (uEfficiency)% if efficiency given
        
        [ax,h1,h2] = plotyy(u,10,uRange,uEfficiency,@myHistogram,@myPlot);
         xlim([uRange(1)-0.05, uRange(end)+0.05]);
         set(get(ax(1), 'YLabel'), 'String', 'freq. of presentation')
         set(get(ax(2), 'YLabel'), 'String', 'efficiency')
         set(ax(1),'YLim', [0 0.4],'YTick',0:.2:.4)       
         xlabel('stimulus intensity')
         box off

    else
        
     % show stimuli density
     histogram(u, 10, 'EdgeColor','none','FaceColor',[.3 .3 .4],'Normalization','probability');
     xlim([uRange(1)-0.05, uRange(end)+0.05])
     ylim([0 0.4])
     xlabel('stimulus intensity')
     ylabel('freq. of presentation')
     box off
    end
     % show type of optimisation
     VBA_title(gca,titleTxt);
       
    if nargin > 6 % if efficiency given
        subplot(4,3,3+idx)
        plot(dEfficiency);
        xlim([1 numel(dEfficiency)])
        ylim([1.1*min(dEfficiency) 0])
        switch idx
            case 2
                xlabel('selected design');
            case 3
                xlabel('trial');
        end

        ylabel('efficiency');
        box off
    end

     % show results if any
     if ~isempty(posterior)
         
        % + observations
        
        subplot(4,3,6+idx)
        % predictions
        opt_plot = options;
        opt_plot.priors = posterior;
        dim_opt = dim;
        dim_opt.p = numel(uRange);
        muy = VBA_getLaplace(uRange',[],@g_psychometric,dim_opt,opt_plot);
        plot(uRange,muy,'r','LineWidth',2);
        % true model
        hold on
        plot(uRange,g_psychometric([],phi,uRange),'Color',[0 .8 0],'LineWidth',2);
        % data
        plot(u,y+jitter(1:numel(y)),'.k');
        
        % options
        ylim([-0.1 1.1])
        xlim([uRange(1)-0.05, uRange(end)+0.05])
        xlabel('stimulus intensity')
        ylabel('prob. of detection')
        hold off
        box off
        if idx == 1
            text(.2,.6,'true model','Color',[0 .8 0]);
            text(.2,.4,'observations','Color','k');
            text(.2,.2,'predicted','Color','r');
        end
     
        % + parameters

        subplot(4,3,9+idx);
        
        % posterior estimates
        plotUncertainTimeSeries(posterior.muPhi,sqrt(diag(posterior.SigmaPhi)),[],gca);
        % true values
        hold on
        plot(phi,'o','MarkerFaceColor',[0 .8 0], 'MarkerEdgeColor',[0 .8 0]);
        % options
        set(gca,'XTickLabel',{'center','slope'})
        ylim([-2 4])
        xlabel('parameter')
        ylabel('posterior estimate')
        hold off
        box off
    end
    
    drawnow
end

    function h = myHistogram (x,y)
        h = histogram(x,y,'EdgeColor','none','FaceColor',[.3 .3 .4],'Normalization','probability');
    end

    function h = myPlot (x,y)
        h = plot(x,y);
        hold on
        [mE, iE] = max (y);
        plot(uRange(iE),mE,'o','MarkerFaceColor',[.3 .3 .3],'MarkerEdgeColor',[.3 .3 .3])
        hold off
    end
end
