function [g,dgdx,dgdP] = g_ttest(x,P,u,in)

% transformations
%%% normal priors on phi_1 
sP = P;
dsPdP = ones(numel(P),1);
%%% sparse laplace priors on phi_2 
[sP(2),dsPdP(2)] = VBA_sparsifyPrior (P(2));


% prediction
g = [sP(1) ; 
     sP(1) + sP(2) ];

% derivatives
dgdP = diag(dsPdP);
dgdP(1,2) = 1*dsPdP(1);
dgdx = [];


% % prediction
% g = [sP(1) - sP(2)/2; 
%      sP(1) + sP(2)/2 ];
% 
% % derivatives
% dgdP = zeros(2);
% dgdP(1,1) = dsPdP(1);
% dgdP(2,1) = -(1/2)*dsPdP(2);
% dgdP(1,2) = dsPdP(1);
% dgdP(2,2) = +(1/2)*dsPdP(2);
% dgdx = [];

end