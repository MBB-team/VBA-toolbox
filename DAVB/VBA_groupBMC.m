function [posterior,out] = VBA_groupBMC(L,options)
% VB algorithm for group-level Bayesian Model Comparison
% function [posterior,out] = VBA_groupBMC(L,options)
% IN:
%   - L: Kxn matrix of log-model evidences (K models; n subjects)
%   - options: a structure containing the following fields:
%       .priors: this variable can contain a field .a, which is the Kx1
%       vector of dummy 'prior counts' of each model
%       .MaxIter: max number of iterations
%       .MinIter: min number of iterations
%       .TolFun: max change in Free Energy

[K,n] = size(L);

%-- fill in options with defaults if needed
if ~exist('options','var')
    options = [];
end
if ~isfield(options,'priors')
    priors = [];
else
    priors = options.priors;
end
if ~isfield(options,'MaxIter')
    options.MaxIter = 32;
end
if ~isfield(options,'MinIter')
    options.MinIter = 1;
end
if ~isfield(options,'TolFun')
    options.TolFun = 1e-2;
end
if ~isfield(options,'DisplayWin')
    options.DisplayWin = 1;
end
if isempty(priors) || ~isfield(priors,'a') || isempty(priors.a)
    priors.a = 1e0*ones(K,1);
    options.priors = priors;
end

%-- initialize posterior and free energy
f0 = priors.a./sum(priors.a);
priors.r = repmat(f0,1,n);
posterior = priors;
F = FE(L,posterior,priors);

%-- enter VB iterative algorithm
stop = 0;
it = 1;

while ~stop
   
    it = it+1;
    
    % update subject attributions
    Elogr = psi(posterior.a) - psi(sum(posterior.a));
    for i=1:n
        tmp = L(:,i)+Elogr;
        g = exp(tmp-max(tmp));
        posterior.r(:,i) = g./sum(g);
    end
    
    % update posterior model counts
    posterior.a = priors.a + posterior.r*ones(n,1);
    
    % calculate Free Energy
    F(it) = FE(L,posterior,priors);
    
    % check stopping criteria
    stop = checkStop(it,F,options);
    
    if options.DisplayWin
        out = wrapUp(L,posterior,priors,F,options);
        options.handles = VBA_displayGroupBMC(posterior,out);
    end
    
end

%-- wrap up VBA output
out = wrapUp(L,posterior,priors,F,options);



function out = wrapUp(L,posterior,priors,F,options)
% wraps up the ou structure for display purposes
[K,n] = size(L);
out.options = options;
out.L = L;
out.F = F;
a0 = sum(posterior.a);
% derive first and second order moments on model frequencies:
out.Ef = posterior.a./a0;
out.Vf = -posterior.a*posterior.a';
out.Vf = out.Vf - diag(diag(out.Vf)) + diag(posterior.a.*(a0-posterior.a));
out.Vf = out.Vf./((a0+1)*a0^2);
% store accuracy and entropy terms of the Free Energy
[F,out.ELJ,out.Sqf,out.Sqm] = FE(L,posterior,priors);
% derive Free Energy under the null:
[out.F0] = FE_null(L);
% Let df(k) = f(k)-max(f(\k)). 
% Then: P(f(k)>f(\k)) = P(df(k)>0)
% Thus, the exceedance probability can be approximated as follows:
% 1-Approximate the random variables f with a Gaussian moment-matching
%   procedure.
% 2-For each k, derive the first- and second- order moments of df(k)
% 3-Calculate P(df(k)>0)
out.ep = zeros(K,1);
c = [1;-1];
for k=1:K
    m = [out.Ef(k);0];
    tmp = out.Ef;
    tmp(k) = -Inf;
    [m(2),im] = max(tmp);
    v = c'*out.Vf([k;im],[k;im])*c;
    mdf = c'*m;
    [out.ep(k)] = VB_PPM(mdf,v,0,0);
end
out.ep = out.ep./sum(out.ep);


function stop = checkStop(it,F,options)
% checks stopping criteria
stop = 0;
if it<options.MinIter
    return
end
dF = F(it) - F(it-1);
if abs(dF)<=options.TolFun || it>=options.MaxIter
    stop = 1;
end


function [F,ELJ,Sqf,Sqm,F0] = FE(L,posterior,priors)
% derives the free energy for the current approximate posterior
[K,n] = size(L);
a0 = sum(posterior.a);
Sqf = sum(gammaln(posterior.a)) - gammaln(a0) + (a0-K)*psi(a0) - sum((posterior.a-1).*psi(posterior.a));
Sqm = 0;
for i=1:n
    Sqm = Sqm - sum(posterior.r(:,i).*log(posterior.r(:,i)));
end
Elogr = psi(posterior.a) - psi(sum(posterior.a));
ELJ = -sum(gammaln(priors.a)) + gammaln(sum(priors.a));
for k=1:K
    ELJ = ELJ + (priors.a(k)-1)*Elogr(k);
    for i=1:n
        ELJ = ELJ + posterior.r(k,i).*(Elogr(k)+L(k,i));
    end
end
F = ELJ + Sqf + Sqm;


function [F0] = FE_null(L)
% derives the free energy of the 'null' (H0: equal model frequencies)
[K,n] = size(L);
F0 = 0;
for i=1:n
    tmp = L(:,i) - max(L(:,i));
    g = exp(tmp)./sum(exp(tmp));
    for k=1:K
        tmp = L(:,i);
        cst = max(tmp);
        tmp = tmp -cst;
        F0 = F0 + g(k).*(log(sum(exp(tmp)))+cst-log(K));
    end
end




