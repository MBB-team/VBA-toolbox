function [py,gridy] = empiricalHist(y)
% numerical derivation of the empirical distribution of y
% function  [py,gridy] = empiricalHist(y)
% IN:
%   - y: nXp data matrix. the histograms are computed along the columns of
%   y, i.e. each column is considered as a vector of samples
% OUT:
%   - py: (n-1)xp matrix of estimated emprical histograms
%   - gridy: (n-1)xp matrix of grid over which the emprical histograms are
%   estimated.
% Note: the following lines of code plots the empirical histogram:
% > [py,gridy] = empiricalHist(y);
% > figure,plot(gridy,py)

[n,p] = size(y);
D = (diag(ones(n,1),0)+diag(ones(n-1,1),1))./2;
D(n,:) = [];
kernel = exp(-(-n/2:n/2).^2./n);
kernel = kernel(:)./sum(kernel);
for i=1:p
    sy = sort(y(:,i));
    [sy,ecdf] = unique(sy);
    dy = (sy(end) - sy(1))./(n-1);
    gridyi = [sy(1):dy:sy(end)]';
    ecdf = interp1(sy,ecdf./n,gridyi);
    py(:,i) = conv(diff(ecdf),kernel,'same');
    gridy(:,i) = D*gridyi;
end
    
