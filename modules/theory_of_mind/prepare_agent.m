function [info] = prepare_agent(model,game,role)
% set basic params for different sorts of agents engaging in dyadic games
% function [info] = prepare_agent(model,game,role)
% info: structure containing the following fields:
%   .f_p = handle of evolution function
%	.g_p = handle of observation function
%   .theta = evolution params 
%   .phi = observation params
%	.x0 = initial conditions
%	.inG = optional inputs to the observation function
%	.inF = optional inputs to the evolution function


switch model
    case 'RB'
        f = @f_Id;
        g = @g_sig;
        theta = [];
        if rand>0.5
            phi = [0;VBA_sigmoid(0.65, 'inverse', true)]; % bias =65%
        else
            phi = [0;VBA_sigmoid(0.35, 'inverse', true)]; % bias =35%
        end
        inF = struct('game',game,'player',role);
        inG = [];
        x0 = 0;
    case 'WSLS'
        f = @f_wslsinGame;
        g = @g_softmax;
        theta = [];
        phi = [1;0]; % inverse (-log) temperature & bias
        inF = struct('game',game,'player',role);
        inG = [];
        x0 = [0;0];
    case 'RL'
        f = @f_RLinGame;
        g = @g_softmax;
        theta = 0; % learning rate = 0.5
        phi = [1;0]; % inverse (-log) temperature & bias
        inF = struct('game',game,'player',role);
        inG = [];
        x0 = [0;0];
    case 'BSL'
        f = @f_BSLinGame;
        g = @g_BSLinGame;
        theta = -2;%[-log(2)]; % (log-) prior volatility
        phi = [-1;0]; % (log-) temperature & bias
        inF.K = 1; % sequence length
        inG = struct('K',inF.K,'game',game,'player',role);
        x0 = zeros(2^(inF.K+1),1);
    case 'HGF'
        f = @f_HGFinGame;
        g = @g_HGFinGame;
        theta = [0;-2;0];%[0;-2;0];
        phi = [-1;0]; % (log-) temperature & bias
        inF.lev2 = 1; % 3rd level (volatility learning)
        inF.kaub = 1.4;
        inF.thub = 1;
        inF.rf = 1;
        inG = struct('game',game,'player',role);
        x0 = [0;0;0;0;0];
    case 'Inf'
        f = @f_Hampton;
        g = @g_Hampton;
        theta = [VBA_sigmoid(0.25,'inverse',true);VBA_sigmoid(0.25,'inverse',true);0]; % (invsigmoid-) weight (PE1), (invsigmoid-) weight (PE2), (log-) opponent's temp
        phi = [-1;0]; % (-log) temperature and bias
        inF = struct('game',game,'player',role);
        inG = struct('game',game,'player',role);
        x0 = 0;
    case '0-ToM'
        K = 0; % depth of k-ToM's recursive beliefs
        [options,dim] = prepare_kToM(K,game,role,0);
        f = @f_kToM;
        g = @g_kToM;
        theta = -2;%-log(2); % (log-) prior volatility
        phi = [-1;0]; % (log-) temperature and bias
        inF = options.inF;
        inG = options.inG;
        x0 = options.priors.muX0;
    case '1-ToM'
        K = 1; % depth of k-ToM's recursive beliefs
        [options,dim] = prepare_kToM(K,game,role,0);
        f = @f_kToM;
        g = @g_kToM;
        theta = -2;%-log(2); % (log-) prior volatility
        phi = [-1;0]; % (log-) temperature and bias
        inF = options.inF;
        inG = options.inG;
        x0 = options.priors.muX0;
    case '2-ToM'
        K = 2; % depth of k-ToM's recursive beliefs
        diluteP = 0; % partial forgetting of opponent's level
        [options,dim] = prepare_kToM(K,game,role,diluteP);
        f = @f_kToM;
        g = @g_kToM;
        theta = -2;%-log(2); % (log-) prior volatility
        phi = [-1;0]; % (log-) temperature and bias
        inF = options.inF;
        inG = options.inG;
        x0 = options.priors.muX0;
    case '3-ToM'
        K = 3; % depth of k-ToM's recursive beliefs
        diluteP = 0; % partial forgetting of opponent's level
        [options,dim] = prepare_kToM(K,game,role,diluteP);
        f = @f_kToM;
        g = @g_kToM;
        theta = -2;%-log(2); % (log-) prior volatility
        phi = [-1;0]; % (log-) temperature and bias
        inF = options.inF;
        inG = options.inG;
        x0 = options.priors.muX0;
    case 'metaToM'
        kToM_level = 1;
        seq_length = 1;
        [theta,phi,inF,inG,x0] = prepare_metaToM(kToM_level,seq_length,game,role);
        f = @f_metaToM;
        g = @g_metaToM;
end

% wrap-up
info.f_p = f;
info.g_p = g;
info.theta = theta; 
info.phi = phi;
info.x0 = x0;
info.inG = inG;
info.inF = inF;
info.payoffTable = game;
