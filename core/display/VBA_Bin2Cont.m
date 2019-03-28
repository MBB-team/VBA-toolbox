function [stacky,stdy,gridg,stdg] = VBA_Bin2Cont(gx,y,maxn)
% transform binary data into continuous intervals (for display purposes)
% function [stacky,stdy,gridg] = VBA_Bin2Cont(gx,y,maxn)
% This function first partitions the data according to percentiles of the
% empirical disitrbution of conditional likelihoods. Then, within each of
% these cells, it calculates the mean and standard deviation of the data.
% IN:
%   - gx: the matrix of conditional likelihoods p(y=1|P,m)
%   - y: the binary data matrix
%   - maxn: the maximum number of cells composing the partition
% OUT:
%   - stacky/stdy: the average (standard deviation, resp.) of data points
%   lying within each partition cell
%   - gridg/stdy: the average (standard deviation, resp.) of conditional
%   likelihoods within each partition cell

y = VBA_vec(y(~isnan(y)));
gx = VBA_vec(gx(~isnan(gx)));

if isempty(gx)
    stacky = [];
    stdy = [];
    gridg = [];
    stdg = [];
    return
end

ny = numel(y);
try;maxn;catch;maxn=min([floor(ny/2),8]);end
ne = min([maxn,ny]);
p = 0:(1/ne):1;
sg = sort(gx);
ind = ceil(p*length(sg)) + [1,zeros(1,size(p,2)-1)];
edges = sg(ind);

stacky = zeros(ne,1);
stdy = zeros(ne,1);
gridg = zeros(ne,1);
stdg = zeros(ne,1);
for i=1:ne
    ind = find(gx>=edges(i)&gx<=edges(i+1));
    stacky(i) = mean(y(ind));
    stdy(i) = std(y(ind));
    gridg(i) = mean(gx(ind));
    stdg(i) = std(gx(ind));
end


