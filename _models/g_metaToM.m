function [gx] = g_metaToM(x,P,u,inG)
% observation function for meta-learner (k-ToM vs sequence)
% function [gx] = g_metaToM(x,P,u,inF)
% The meta-learner bases her decision (a=1 or a=0) upon her prediction of
% her opponent's next move, given the game payoff table. The specific
% difficulty of the meta-learner is that she does not know with absolute
% certainty whether she faces an intentional agent (k-ToM) or an inanimate
% patter (sequence). In fact, she holds a belief about this, which is
% quantified in terms of the probability Pi that she is facing an
% intentional agent (the probability that the agent faces an inanimate
% pattern is 1-Pi).
% IN:
%   - x: hidden states (see indexing in inF.indlev)
%   - P: observation params: phi(1) = log-temperature and phi(2) = bias
%   - u: sequence of past actions:
%       u(1)= opponent's last move
%       u(2)= learner's last move
%       u(3:K+2) = sequence of K past opponent's moves
%   - inG: input structure (see prepare_metaToM.m)
% OUT:
%   - gx: updated hidden states
% [see RecToMfunction.m and f_BSL.m]


% 1- get k-ToM probabilistic decision P(a=1|k-ToM)
xktom = x(inG.ktom.indx);
inG_ktom = inG.ktom.inG;
gx_ktom = ObsRecGen(xktom,P,u,inG_ktom);

% 2- get BSL probabilistic decision P(a=1|seq)
xseq = x(inG.seq.indx);
inG_seq = inG.seq.inG;
gx_bsl = g_BSLinGame(xseq,P,u,inG_seq);


% 3- derive meta-ToM probabilistic decision P(a=1)
Pi = VBA_sigmoid(x(inG.meta.indx)); % prior P(agent=kToM)
gx = Pi*gx_ktom + (1-Pi)*gx_bsl;







