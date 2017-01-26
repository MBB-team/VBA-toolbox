function [options,dim] = prepare_kToM(level,payoffTable,role)
% derives input structures for k-ToM learning models
% function [inF,inG] = prepare_kToM(level,payoffTable)
% IN:
%   - level: depth of k-ToM's recursive beliefs (k)
%   - payoffTable: 2x2x2 payoff table (last dim = player role)
%   - role: player's role (cf. payoff table)
% OUT:
%   - options: structure containing the following fields:
%       .inF/inG: optional input structures for f_kToM.m and g_kToM.m
%       .priors: here, for initial conditions (X0)
%   - dim: model dimensions

inF.player=role; % player role (here: 1=seeker, 2=hider) 
inF.lev=level;
inF.game=payoffTable;
inF.fobs=@ObsRecGen;
inF.indParev=1; % Nb para evol
inF.indParobs=2; %Nb para obs
inF.dummyPar=[1;0;0];
NtotPar=inF.indParev +inF.indParobs;
inG=inF;
inG.indlev= defIndlev(level, NtotPar);
inG.npara=NtotPar;
options.inF = inF;
options.inG = inG;

dim.n = sizeXrec(level,NtotPar);
dim.n_phi = 2;
dim.n_theta = 1;

options.priors.muX0 = f_kToM(zeros(dim.n,1),zeros(2,1),[],inF);
options.priors.SigmaX0 = zeros(dim.n);