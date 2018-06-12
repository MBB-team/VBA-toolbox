function  [ gx,dgdx,dgdP ] = g_softmax(x,P,u,in )
% softmax decision rule for Q-learning (2-armed bandit task)
% function  [ gx,dgdx,dgdP ] = g_softmax(x,P,u,in )
% IN:
%   - x : Q-values
%   - P : inverse (log-) temperature and bias
%   - u : [useless]
%   - in : [useless]
% OUT:
%   - gx : P(a=1|x)

% inverse temperature
% -------------------------------------------------------------------------
beta = exp(P(1)); % exp: [-Inf,Inf] -> [0 Inf]

dQ = (x(1)-x(2));
if length(P)>1
    gx = VBA_sigmoid( beta*dQ + P(2));
else
    gx = VBA_sigmoid( beta*dQ );
end
dgdx = zeros(size(x,1),1);
dgdx(1) = beta*gx*(1-gx);
dgdx(2) = -beta*gx*(1-gx);
if length(P)>1
    dgdP = [beta*dQ*gx*(1-gx),gx*(1-gx)];
else
    dgdP = [beta*dQ*gx*(1-gx)];
end
