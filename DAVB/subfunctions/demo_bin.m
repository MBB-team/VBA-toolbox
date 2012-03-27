% demo for binary data classification
%------------------------------------------------------------
% Copyright (C) 2012 Jean Daunizeau / License GNU GPL v2
%------------------------------------------------------------

close all
clear all

dim.p = 5e2;
dim.n_t = 1;
dim.n_phi = 8;
dim.n_theta = 0;
dim.n = 0;

phi = randn(dim.n_phi,1);
g_fname = @g_classif;


options.inG.X = randn(dim.n_phi-1,dim.p);
options.binomial = 1;
options.priors.muPhi = zeros(dim.n_phi,1);
options.priors.SigmaPhi = 1e0*eye(dim.n_phi);



[y] = simulateNLSS(dim.n_t,[],g_fname,[],phi,[],[],[],options,[]);

% options.checkGrads = 1;

options.isYout = randn(size(y))>2;
[posterior,out] = VBA_NLStateSpaceModel(y,[],[],g_fname,dim,options);

displayResults(posterior,out,y,[],[],[],phi,[],[])



nmcmc = 2e4;
q = zeros(nmcmc,1);
mu = options.priors.muPhi;
sS = getISqrtMat(options.priors.SigmaPhi,0);
et0 = clock;
fprintf(1,'MCMC estimate of the model evidence...')
fprintf(1,'%6.2f %%',0)
for i=1:nmcmc
    Pi = mu + sS*randn(dim.n_phi,1);
    gi = feval(g_fname,[],Pi,[],options.inG);
    tmp = [gi.^y,(1-gi).^(1-y)];
    q(i) = prod(tmp(:));
    if mod(100*i/nmcmc,10) <1
        fprintf(1,repmat('\b',1,8))
        fprintf(1,'%6.2f %%',floor(100*i/nmcmc))
    end
end
fprintf(1,repmat('\b',1,8))
fprintf(1,[' OK (took ',num2str(etime(clock,et0)),' seconds).'])
fprintf(1,'\n')
lpy0 = log(mean(q))
% 
% X = [options.inG.X',ones(size(options.inG.X,2),1)];
% tmp = X*options.priors.muPhi;
% tmp = tmp./(sqrt(1+0.368.*diag(X*options.priors.SigmaPhi*X')));
% Es = 1./(1+exp(-tmp));
% lpy = sum(log(1-y+Es.*(-1).^(1-y)))

