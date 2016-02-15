function [p_BMA] = VBA_BMA(p0,F0)
% performs Bayesian Model Averaging (BMA)
% function [p,o] = VBA_BMA(posterior,F)
% IN:
%   - p0: a Kx1 cell-array of VBA posterior structures, which are
%   conditional onto specific generative models (where K is the number of
%   models)
%   - F0: a Kx1 vector of log-model evidences
% OUT:
%   - p_BMA: the resulting posterior structure, with the first two moments of
%   the marginal probability density functions

K = length(p0); % # models
ps = softmax(F0);

% observation parameters
mus = cell(K,1);
Qs = cell(K,1);
for k=1:K
    mus{k} = p0{k}.muPhi;
    Qs{k} = p0{k}.SigmaPhi;
end
[p_BMA.muPhi,p_BMA.SigmaPhi] = get2moments(mus,Qs,ps);

% evolution parameters
mus = cell(K,1);
Qs = cell(K,1);
for k=1:K
    mus{k} = p0{k}.muTheta;
    Qs{k} = p0{k}.SigmaTheta;
end
[p_BMA.muTheta,p_BMA.SigmaTheta] = get2moments(mus,Qs,ps);

% initial conditions
mus = cell(K,1);
Qs = cell(K,1);
for k=1:K
    mus{k} = p0{k}.muX0;
    Qs{k} = p0{k}.SigmaX0;
end
[p_BMA.muX0,p_BMA.SigmaX0] = get2moments(mus,Qs,ps);

% hidden states
T = size(p0{1}.muX,2);
for t=1:T
    mus = cell(K,1);
    Qs = cell(K,1);
    for k=1:K
        mus{k} = p0{k}.muX(:,t);
        Qs{k} = p0{k}.SigmaX.current{t};
    end
    [p_BMA.muX(:,t),p_BMA.SigmaX.current{t}] = get2moments(mus,Qs,ps);
end

% data precision
mus = cell(K,1);
Qs = cell(K,1);
for k=1:K
    mus{k} = p0{k}.a_sigma/p0{k}.b_sigma;
    Qs{k} = p0{k}.a_sigma/p0{k}.b_sigma^2;
end
[m,v] = get2moments(mus,Qs,ps);
p_BMA.b_sigma = m/v;
p_BMA.a_sigma = m*p_BMA.b_sigma;

% hidden state precision
mus = cell(K,1);
Qs = cell(K,1);
id = zeros(K,1);
for k=1:K
    mus{k} = p0{k}.a_alpha/p0{k}.b_alpha;
    Qs{k} = p0{k}.a_alpha/p0{k}.b_alpha^2;
    if isinf(p0{k}.a_alpha) && p0{k}.b_alpha==0
        id(k) = 1;
    end
end
if isequal(sum(id),K) % all deterministic systems
    p_BMA.b_alpha = Inf;
    p_BMA.a_alpha = 0;
elseif isequal(sum(id),0) % all stochastic systems
    [m,v] = get2moments(mus,Qs,ps);
    p_BMA.b_alpha = m/v;
    p_BMA.a_alpha = m*p_BMA.b_sigma;
else
    disp('VBA_MBA: Warning: mixture of deterministic and stochastic models!')
    p_BMA.b_alpha = Inf;
    p_BMA.a_alpha = 0;
end






function [m,V] = get2moments(mus,Qs,ps)
V = zeros(size(Qs{1}));
m = zeros(size(mus{1}));
K = length(ps);
for k=1:K
    m = m + ps(k).*mus{k};
end
for k=1:K
    tmp = mus{k} - m;
    V = V + ps(k).*(tmp*tmp' + Qs{k});
end