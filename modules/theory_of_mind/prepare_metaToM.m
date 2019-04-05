function [theta,phi,inF,inG,x0] = prepare_metaToM(kToM_level,seq_length,game,role)
% set basic params for meta-ToM learner engaging in dyadic games
% function [f,g,theta,phi,inF,inG,x0] = prepare_metaToM(kToM_level,seq_length,game,role)
% 

theta = [-2;-2]; % prior log-volatility (k-ToM & BSL)
phi = [-1;0]; % (log-) temperature and bias

inF.meta.indx = 1; % (invsigmoid-) P(opponent=k-ToM)
inF.meta.indP = [];
inF.meta.diluteP = 0;
inG.meta.indx = inF.meta.indx;
x0(1) = 0; % P(opponent=k-ToM)=1/2

% prepare kToM
[options,dim] = prepare_kToM(kToM_level,game,role,0);
inF.ktom.inF = options.inF;
inG.ktom.inG = options.inG;
inF.ktom.indx = 2:dim.n+1; % k-ToM states
inF.ktom.indP = length(inF.meta.indP)+1;
inG.ktom.indx = inF.ktom.indx;
x0(inF.ktom.indx) = options.priors.muX0;

% prepare BSL (sequence learner)
inF.seq.inF.K = seq_length; % sequence length
inG.seq.inG = struct('K',seq_length,'game',game,'player',role);
inF.seq.indx = inF.ktom.indx(end)+1:inF.ktom.indx(end)+2^(seq_length+1);
inF.seq.indP = inF.ktom.indP+1;
inG.seq.indx = inF.seq.indx;
x0(inF.seq.indx) = zeros(2^(seq_length+1),1);
x0 = VBA_vec(x0);



