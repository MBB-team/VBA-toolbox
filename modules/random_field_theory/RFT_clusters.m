function [clusters,imax] = RFT_clusters(X,xc,verbose)
% identifies the clusters on a thresholded RF defined on a regular lattice
% [out] = RFT_clusters(X,xc,verbose)
% This function finds the connected clusters of points Pi on a regular
% lattice, such that their value X(Pi) exceeds the threshold xc.
% IN:
%   - X: the nX1 vector or discretized RF on the regular lattice
%   - xc: the cluster-inducing threshold
%   - verbose: verbose flag
% OUT:
%   - clusters: a ncX1 cell array containing the indices of the lattice
%   vertices belonging to each cluster, where nc is the number of connected
%   clusters induced by xc.
%   - imax: a ncX1 array of indices of clusters' maxima


try,verbose;catch,verbose=0;end
X = VBA_vec(X);
% 1- find upcrossing clusters induced by xc
B = [0;VBA_vec(X>xc);0];
n = length(X);
in = VBA_vec(find(B==1));
if isempty(in)
    clusters = [];
    nc = 0;
else
    dB = [diff(B);0];
    tmp = cumsum(dB(in-1));
    C = zeros(size(B));
    C(in) = tmp;
    nc = max(tmp);
    clusters = cell(nc,1);
    for i=1:nc
        clusters{i} = find(C==i)-1;
    end
end
% 2-find clusters' maxima
imax = [];
for i=1:nc
    [xm,im] = max(X(clusters{i}));
    imax(i) = clusters{i}(im(1));
end

if ~verbose
    return
end
col = getColors(nc);
hf = figure('color',[1 1 1]);
ha = axes('parent',hf,'nextplot','add');
plot(ha,X,'k','marker','.')
plot(ha,[1,n],[xc,xc],'k--')
for i=1:nc
    plot(ha,clusters{i},X(clusters{i}),'.','color',col(i,:))
    plot(ha,imax(i),X(imax(i)),'o','color',col(i,:))
end

