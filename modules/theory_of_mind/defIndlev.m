function indlev = defIndlev(level,NtotPar)
% indexing of k-ToM's hidden states
% function findlev = defIndlev(level,NtotPar)
% IN:
%   - level: k-ToM's sophistication level
%   - NtotPar: total nb of evol/obs params
% OUT:
%   - indlev: structure of hidden-states indices

if level==0 % 0-ToM
    indlev = [];
else % k-ToM with k>0
    curind = level-1;
    for i=1:level % loop over k-ToM admissible opponents
        indlev(i).k = i-1; % opponent's sophistication
        sizeXi = sizeXrec(i-1,NtotPar);
        indlev(i).X = curind+(1:sizeXi); % opponent's hidden states
        curind = curind+sizeXi;
        indlev(i).f = curind+1; % x(theta) = log odds of P(o=1|k'=i-1)
        curind = curind+1;
        indlev(i).df = curind+(1:NtotPar); % dx/dtheta
        curind = curind+NtotPar;
        indlev(i).Par = curind+(1:2*NtotPar); % E[theta|k'=i-1] and V[theta|k'=i-1]
        curind = curind+2*NtotPar;
    end
end
