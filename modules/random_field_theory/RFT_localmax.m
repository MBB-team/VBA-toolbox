function imax = RFT_localmax(X)
% find local maxima on the field
% function imax = RFT_localmax(X)
% This function searches for local maxima on the field, whose values X are
% assumed to be sampled on a regular lattice.
% IN:
%   - X: the Lx1 sampled RF
% OUT:
%   - imax: set of locam maxima indices

L = length(X);

% D1 = diag(ones(L,1),0)-diag(ones(L-1,1),1);
% D2 = diag(ones(L,1),0)-diag(ones(L-1,1),-1);
% imax = vec(find(D1*X>0&D2*X>0));

X1 = [-Inf;X];
X2 = [X;-Inf];
imax = intersect(find(diff(X1)>0),find(diff(X2)<0));
imax = VBA_vec(imax);

if ~isequal(imax(1),1) && X(1)>X(2)
    imax = [1;imax];
end
if ~isequal(imax(end),L) && X(L)>X(L-1)
    imax = [imax;L];
end

% if ~isequal(imax2(1),1) && X(1)>X(2)
%     imax2 = [1;imax2];
% end
% if ~isequal(imax2(end),L) && X(L)>X(L-1)
%     imax2 = [imax2;L];
% end
% isequal(imax,imax2)
% figure
% plot(X)
% hold on
% plot(imax,X(imax),'o')
% plot(imax2,X(imax2),'+')
% pause