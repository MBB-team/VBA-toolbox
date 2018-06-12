function [fx] = f_metaToM(x,P,u,inF)
% evolution function for meta-learning (k-ToM vs sequence)
% function [fx] = f_metaToM(x,P,u,inF)
% This function update the belief of an agent who has to identify whether
% he faces an intentional agent (k-ToM) or an inanimate patter (sequence).
% Learning here derives from VB-Laplace update rules of the agent's
% posterior belief. practically speaking, this function (i) wraps modified
% evolution functions of both k-ToM learning and sequence learning, and
% (ii) updates the probability Pi that the agent is facing an intentional
% agent given a new observation u (the probability that the agent faces an
% inanimate pattern is 1-Pi).
% IN:
%   - x: hidden states (see indexing in inF.indlev)
%   - P: evolution params:
%   P(1)= agent's prior opponent's (log-) volatility on opponent's params
%   P(2)= agent's (invsigmoid-) dilution coefficient
%   - u: sequence of past actions:
%       u(1)= opponent's last move
%       u(2)= learner's last move
%       u(3:K+2) = sequence of K past opponent's moves
%   - inF: input structure (see prepare_metaToM.m)
% OUT:
%   - fx: updated hidden states
% [see RecToMfunction.m and f_BSL.m]

ot = u(1); % opponent's last move
fx = NaN(size(x)); % initialize updated states

% 1- update P(agent=kToM)
Pi0 = VBA_sigmoid(x(inF.meta.indx)); % prior P(agent=1)
% partial forgetting of prior belief on opponent's type?
if inF.meta.diluteP
    dc = VBA_sigmoid(P(inF.meta.indP)); % dilution coefficient
    Pi0 = (1-dc).*Pi0 + dc./2;
end

% Get k-ToM's likelihood, ie P(o|k-ToM).
xktom = x(inF.ktom.indx);
inFktom = inF.ktom.inF;
ntotPar = inFktom.indParev+inFktom.indParobs; % total number of params
level = inFktom.lev; % depth of k-ToM's recursive beliefs (k=level)
indlev = defIndlev(level,ntotPar); % states indexing
if level==0 % 0-ToM [should be useless]
    mx = xktom(1); % E[log-odds of P(o=1)]
    Vx = exp(xktom(2)); % V[log-odds of P(o=1)]
    Els1 = VBA_Elogsig(mx,Vx);
    Els0 = VBA_Elogsig(-mx,Vx);
    ELL = ot.*Els1 + (1-ot).*Els0;
    h_kToM = exp(ELL); % P(o|k-ToM)
else
    Pk = VBA_sigmoid(x(1:(level-1))); % P(k'), with k'=0,...,k-1
    Pk = [Pk;max(0,1-sum(Pk))]; % insert last P(k'=k-1)  
    f = zeros(level,1); % E[x(theta)]
    Vx = zeros(level,1); % V[x(theta)]
    for j=1:level % loop over possible opponent's levels  (k'=j-1)
        f(j) = xktom(indlev(j).f); % E[x(theta)|k'=j-1]
        df = xktom(indlev(j).df); % d[x(theta)]/dtheta for k'=j-1
        Sig = exp(xktom(indlev(j).Par(2:2:2*ntotPar))); % V[theta|k'=j-1]
        Vx(j) = sum(Sig.*df.^2); % V[x(theta)|k'=j-1]
    end
    Els1 = VBA_Elogsig(f,Vx);
    Els0 = VBA_Elogsig(-f,Vx);
    ELL = ot.*Els1 + (1-ot).*Els0;
    h_kToM = exp(Pk'*ELL); % P(o|k-ToM)
end

% Get sequence's likelihood, ie P(o|seq).
xseq = x(inF.seq.indx);
inFseq = inF.seq.inF;
useq = u; % remove agent's previous move
useq(2) = [];
K = inFseq.K; % sequence depth
yb = useq(2:K+1); % previous outcomes
if VBA_isWeird (yb)
    h_seq = 1/2;
else
    if K >0
        indSeq = bin2dec(num2str(yb'))+1; % index of sequence of previous outcomes
    else
        indSeq = 1;
    end
    m = xseq(indSeq);
    v = exp(xseq((2^K)+indSeq));
    Els1 = VBA_Elogsig(m,v);
    Els0 = VBA_Elogsig(-m,v);
    ELL = ot.*Els1 + (1-ot).*Els0;
    h_seq = exp(ELL); % P(o|seq)
end

% VB update of P(agent=kToM)
Pi = Pi0.*h_kToM./(Pi0.*h_kToM+(1-Pi0).*h_seq);
fx(inF.meta.indx) = VBA_sigmoid(Pi, 'inverse', true);


% 2- update k-ToM belief
inFktom.metaweight = Pi;
Par_ktom = P(inF.ktom.indP);
uktom = u(1:2);
fx(inF.ktom.indx) = RecToMfunction(xktom,Par_ktom,uktom,inFktom);

% 3- update seq learner belief
inFseq.metaweight = 1-Pi;
Par_seq = P(inF.seq.indP);
fx(inF.seq.indx) = f_BSL(xseq,Par_seq,useq,inFseq);





