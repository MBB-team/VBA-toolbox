function Ltilde = VBA_BTWgroup(L)
% derives subjects' evidence for same/different model across conditions
% function Ltilde = VBA_BTWgroup(L)
% IN:
%   - L: Kxnx2 array of log-model evidences (K models; n subjects; 2
%   conditions).
% OUT:
%   - Ltilde: 2xn array of log-evidence of c0 (same model for both
%   conditions) and c1 (distinct models across conditions).

[K,n,c] = size(L);
if c~=2
    disp('Error: log-model evidences array should be Kxnx2!')
    Ltilde = [];
    return
end
Ltilde = zeros(2,n);
for i=1:n
    LL1 = L(:,i,1);
    LL2 = L(:,i,2);
    mL1 = max(LL1);
    mL2 = max(LL2);
    p1 = exp(LL1-mL1);
    p2 = exp(LL2-mL2);
    % evidence for c0
    e0 = sum(p1.*p2);
    Ltilde(1,i) = mL1 + mL2 + log(eps+e0./K);
    % evidence for c1:
    e1 = sum(p1);
    e2 = sum(p2);
    Ltilde(2,i) = mL1 + mL2 + log(eps+(e1.*e2-e0)./(K*(K-1)));
end