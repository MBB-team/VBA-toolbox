function [priorsMu, priorsSigma,labeler] = priorUglyfier(varargin)



n = numel(varargin);
assert(mod(n,3)==0,'***Parameters should be: ''name'',mu,sigma,...');

priorsMu = [];
priorsSigma = [];
labeler = struct ;

cpt = 0;
for i=1:3:n
    
    name = varargin{i};
    mu = varargin{i+1};
    sigma = varargin{i+2};
    
    
    nPar = numel(mu);
    par_idx = cpt+(1:nPar);
    cpt = cpt + nPar;

    labeler.(name)=par_idx ;
    
    priorsMu(par_idx) = mu;
    
    priorsSigma(par_idx,par_idx) = diag(sigma);

    
    
end
priorsMu = priorsMu(:);
    