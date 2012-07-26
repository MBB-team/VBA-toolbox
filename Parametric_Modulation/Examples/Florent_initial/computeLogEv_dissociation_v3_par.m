function [modev,finite_ind] = computeLogEv_dissociation_v3_par(subList, type, model)
% FOR PARALLEL COMPUTING
% TASK: IMPLICIT
% MODEL: HYPERBOLIQUE & LINEAR
% v3: more models allowed.
% Pipline to compute the Log-evidence of each accumulation model for each
% subject, using the DAVB toolbox
%
% Usage: [modev] = computeLogEv_fn(subList, method, type)
%   IN
%   - subList: list of subject, as call by the extraction function
%   - type: 'FindOnOff' to use extractdata_FindOnOff_fn.m
%           'Markers' to use extractdata_markers_fn.m
%   - model: 'linear', 'interaction' or 'hyperbolic'
%
%   OUT
%   - modev, a nSub*nModel x 9 matrix. Colums are:
%       1: sub #
%       2: model #
%       3: modulation of K by benefits (1: enabled, 0: disabled)
%       4: modulation of Je by benefits (1: enabled, 0: disabled)
%       5: modulation of Jr by benefits (1: enabled, 0: disabled)
%       6: modulation of K by cost (1: enabled, 0: disabled)
%       7: modulation of Je by cost (1: enabled, 0: disabled)
%       8: modulation of Jr by cost (1: enabled, 0: disabled)
%       9: modulation of K by display (1: enabled, 0: disabled)
%       10: modulation of Je by display (1: enabled, 0: disabled)
%       11: modulation of Jr by display (1: enabled, 0: disabled)
%
%Lou Safra 2011-12

% --- SET PATHS ---
% =================

addpath /home/florent/DAVB
addpath /home/florent/DAVB/subfunctions
addpath /home/florent/ACCUMULATEUR/these/manip/ACCU-EXPLICITE-DISSOCIATION/scripts/script-analysis/BMS_lou/BMS_dissociation

% --- PARAMETERS ---
% ==================

% model of observation
switch model
    case 'linear'
        g_fname = @g_RT_dissociation_linear;
        
    case 'hyperbolic'
        g_fname = @g_RT_dissociation_hyperbolic;
end

dim.n_phi   = 12;
dim.n       = 0;
dim.n_theta = 0;

inG.K.mean  = 1;
inG.K.p1    = 2;
inG.K.p2    = 3;
inG.K.p3    = 4;
inG.Je.mean = 5;
inG.Je.p1   = 6;
inG.Je.p2   = 7;
inG.Je.p3   = 8;
inG.Jr.mean = 9;
inG.Jr.p1   = 10;
inG.Jr.p2   = 11;
inG.Jr.p3   = 12;

inG.reward  = zscore([1 1 1 1 2 2 2 2])'+1;
inG.display = zscore([1 1 2 2 1 1 2 2])'+1;
inG.effort  = zscore([1 2 1 2 1 2 1 2])'+1;

% DEFINE BAD MODEL A PRIORI

Kall=[2 6 10 3 7 11 4 8 12];

% StupidIModel 2 6 10
StupidIModel=[0 0 0
              0 0 1
              0 1 0
              1 0 0];

% StupidFModel 3 7 11
StupidFModel=[0 0 0
              0 0 1];

% StupidDModel 4 8 12
StupidDModel=[0 0 0
              0 1 0];

modev = {};

parfor iSub = 1:length(subList)
    
    modev{iSub} = [];
    
    % Set parameter priors
    priors = {}; 
    priors.a_sigma                  = 1e0;                   % Jeffrey's prior
    priors.b_sigma                  = 1e0;                     % Jeffrey's prior
    priors.muPhi                    = 1e-1*ones(dim.n_phi,1);% prior mean on observation params
    priors.SigmaPhi                 = 1e2*eye(dim.n_phi);    % prior covariance on observation params
    
    % options.checkGrads    = 1;
    options = {};
    options.priors          = priors;       % include priors in options structure
    options.inG             = inG;          % input structure (grid)
    options.GnFigs          = 0;            % disable annoying figures
    options.verbose         = 0;
    options.DisplayWin      = 0;
    
    fprintf('\n SUBJECT %d', iSub)
    
    % Get Subject data
    y = getData(subList(iSub), type);

    % --- model  ---
    nmodel=0;
    for i2=[0 1]
        for i6=[0 1]
            for i10=[0 1]
                for i3=[0 1]
                    for i7=[0 1]
                        for i11=[0 1]
                            for i4=[0 1]
                                for i8=[0 1]
                                    for i12=[0 1]
                                        
                                        nmodel = nmodel+1;
                                        Kind=Kall([i6 i10 i3 i7 i11 i4 i8 i12]==0);
                                        
                                         %if isNotStupid
                                        priors.muPhi = 1e-1*ones(dim.n_phi,1);         % prior mean on observation params
                                        priors.SigmaPhi = 1e2*eye(dim.n_phi);          % prior covariance on observation params
                                        
                                        priors.muPhi(Kind) = 0; % parameter set to 0
                                        priors.SigmaPhi(Kind,Kind) = 0;
                                        
                                        
                                        options.priors = priors;
                                        try
                                            doModel = 1;
                                            
                                            switch model
                                                case 'linear'
                                                    for i = 1:size(StupidIModel,1)
                                                        if StupidIModel(i,:)==[i2 i6 i10]
                                                            doModel = 0;
                                                        end
                                                    end
                                                    
                                                    for i=1:size(StupidFModel,1)
                                                        if StupidFModel(i,:)==[i3 i7 i11]
                                                            doModel = 0;
                                                        end
                                                    end
                                                    
                                                    for i=1:size(StupidDModel,1)
                                                        if StupidDModel(i,:)==[i4 i8 i12]
                                                            doModel = 0;
                                                        end
                                                    end
                                                case 'hyperbolic'
                                                    for i = 1:size(StupidIModel,1)
                                                        if StupidIModel(i,:)==[i2 i6 i10]
                                                            doModel = 0;
                                                        end
                                                    end
                                                    
                                                    for i=1:size(StupidFModel,1)
                                                        if StupidFModel(i,:)==[i3 i7 i11]
                                                            doModel = 0;
                                                        end
                                                    end
                                                    
                                                    for i=1:size(StupidDModel,1)
                                                        if StupidDModel(i,:)==[i4 i8 i12]
                                                            doModel = 0;
                                                        end
                                                    end
                                                    
                                                    if i2==0 & i4==1
                                                        doModel =0;
                                                    end
                                                    
                                                    
                                                    if i6==0 & i8==1
                                                        doModel =0;
                                                    end
                                                    
                                                    if i10==0 & i12==1
                                                        doModel =0;
                                                        
                                                    end                                                    
                                            end                                            
                                            
                                            if doModel==1
                                                [p1,o1] = VBA_NLStateSpaceModel(y,[],[],g_fname,dim,options);
                                                
                                                modev{iSub}=[modev{iSub}; iSub nmodel i2 i6 i10 i3 i7 i11 i4 i8 i12 o1.F];
                                                
                                            else
                                                
                                                modev{iSub}=[modev{iSub}; iSub nmodel i2 i6 i10 i3 i7 i11 i4 i8 i12 -Inf];
                                                
                                            end
                                            
                                        catch
                                            
                                            modev{iSub}=[modev{iSub}; iSub nmodel i2 i6 i10 i3 i7 i11 i4 i8 i12 NaN];
                                            
                                            
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end
assignin('base', 'modev', modev)
% sort model after the parfor loop
tmp = [];
for iSub = 1:length(subList)
    tmp = [tmp; modev{iSub}];
end

modev = tmp;

for imodel = unique(modev(:,2))'
    sumev(imodel) = sum(modev(modev(:,2)==imodel, 12));
end
figure(1)
bar(sort(sumev(isfinite(sumev))));
[~, ind] = sort(sumev);
finite_ind = ind(isfinite(sort(sumev)));
set(gca,'XTick',1:length(finite_ind),'XTickLabel',{num2str(finite_ind')})

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function y = getData(iSub, type)
% Compute the subject mean for each RW x FC x DISP condition, using the function
% set by type
%
% Usage: y = getData(iSub, type)
%   - iSub: subID called by the extraction function
%   - type: 'FindOnOff' to use extractdata_FindOnOff_fn.m
%           'Markers' to use extractdata_markers_fn.m
%   - y: the output column vector has effort first and rest last, and the
%       condition in the order :
%       RW FC DISP
%       1  1   1
%       1  1   2
%       1  2   1
%       2  2   2
%       2  1   1
%       2  1   2

switch type
    case 'FindOnOff'
        
        opt.datadir = '/home/florent/ACCUMULATEUR/these/manip/ACCU-EXPLICITE-DISSOCIATION/data-subject/EXPE2B/';
        opt.FilePrefix = 'MotAccSubExpe2b_';
        opt.subIDlist = iSub;
        opt.def_blocks = 'default'; % set number or 'default'
        opt.def_hands = ['d' 'g'];
        opt.def_trials = 'default'; % set number or 'default'
        opt.ON2OFFDeadtime = 0;
        opt.OFF2ONDeadtime = 0;
        opt.thd = 0.5;
        opt.SR = 20;
        opt.sd_threshold = 1;
        opt.method = 'full+coin';  % 'full', 'min', 'min+end', 'min+start' , 'min+start+coin', 'full+coin'
        
        [EffortDur, RestDur, ~, ~] = extractdata_FindOnOff_fn(opt);
      
    case 'Markers'
        
        opt.H.Fs            = 1250;
        opt.H.nSamplesPre   = 1250;
        opt.datadir         ='/mnt/data/ACC_MEG/data/';
        opt.FilePrefix      = 'BLINKCor_';
        opt.subIDlist       = iSub;
        opt.progressbar     = 0; % 1 to display progress bar, 0 not to display it
        
        [EffortDur, RestDur] = extractdata_markers_fn(opt);
end

i=1;
for iRW = unique(EffortDur(:, 5))'
    for iDISP = unique(EffortDur(:,7))'
        for iFC = unique(EffortDur(:,6))'
            y(i+8) = mean(RestDur(RestDur(:,5) == iRW &...
                RestDur(:,6) == iFC &...
                RestDur(:,7) == iDISP, 8));
            
            y(i) = mean(EffortDur(EffortDur(:,5) == iRW &...
                EffortDur(:,6) == iFC &...
                EffortDur(:,7) == iDISP, 8));
            i=i+1;
            
        end
    end    
end
    y=y';
end