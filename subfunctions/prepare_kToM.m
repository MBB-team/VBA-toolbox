function [options,dim] = prepare_kToM(K,payoffTable,role,diluteP)
% specify input structures for k-ToM learning models
% function [inF,inG] = prepare_kToM(level,payoffTable)
% Note: k-ToM models can be extended with partial "forgetting" of the
% opponent's sophisitication level [see RecToMfunction.m].
% IN:
%   - K: depth of k-ToM's recursive beliefs
%   - payoffTable: 2x2x2 payoff table (last dim = player role)
%   - role: player's role (cf. payoff table)
%   - diluteP: flag for partial "forgetting" effect on opponent's level (0:
%   no forgetting, 1: partial forgetting)
% OUT:
%   - options: structure containing the following fields:
%       .inF/inG: optional input structures for f_kToM.m and g_kToM.m
%       .priors: here, for initial conditions (X0)
%   - dim: model dimensions

try,diluteP;catch,diluteP=0;end

inF.player = role; % player role (here: 1=seeker, 2=hider) 
inF.lev = K; % player's sophistication level
inF.game = payoffTable; % game payoff table
inF.fobs = @ObsRecGen; % opponent's observation function 

%--- change the following 3 lines if 0-ToM evol/obs functions change! ---%
inF.diluteP = diluteP; % partial forgetting of opponent's level
if inF.diluteP==0
    inF.indParev = 1; % nb of evol params
    inF.dummyPar = [1;0;0]; % which paramaters are time-dependent
elseif inF.diluteP==1
    inF.indParev = 2; % nb of evol params
    inF.dummyPar = [1;0;0;0]; % which paramaters are time-dependent
end
inF.indParobs = 2; % nb of obs params

NtotPar = inF.indParev +inF.indParobs; % total nb of params
inG = inF;
inG.indlev = defIndlev(K,NtotPar); % states' indexing
inG.npara = NtotPar;
options.inF = inF;
options.inG = inG;

dim.n = sizeXrec(K,NtotPar);
dim.n_phi = inF.indParobs;
dim.n_theta = inF.indParev;

options.priors.muX0 = f_kToM(zeros(dim.n,1),zeros(2,1),[],inF);
% specify initial belief Re: opponent's params
if K > 0 % for k-ToM with k>0
    for j=1:K
        % prior (log-) volatility = -1 (opponent = efficient learner)
        indVol = inG.indlev(j).Par(1); % index of hidden-state encoding E[invsigm volatility]
        options.priors.muX0(indVol) = 0;%-1;
        % "forgetting" of opponent's sophistication
        if inF.diluteP==1
            indDil = inG.indlev(j).Par(3); % index of hidden-state encoding E[invsigm dilut. coef.]
            options.priors.muX0(indDil) = -2;
        end
        % behavioural temperature = -1 (opponent = exploitative decider)
        indTemp = inG.indlev(j).Par(2*inF.indParev+1); % index of hidden-state encoding E[log temp]
        options.priors.muX0(indTemp) = 0;%-1;
    end
end
% options.priors.SigmaX0 = zeros(dim.n);


