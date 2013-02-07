function [py,gridy] = empiricalHist(y,pr)
% numerical derivation of the empirical distribution of y
% function  [py,gridy] = empiricalHist(y)
% IN:
%   - y: nXp data matrix. The histograms are derived along the columns of
%   y, i.e. each column is considered as a vector of samples
%   - pr: the precision of the smoothing (gaussian) kernel {1}.
% OUT:
%   - py: (n-1)xp matrix of estimated emprical histograms
%   - gridy: (n-1)xp matrix of grid over which the emprical histograms are
%   estimated.
% Note: the following lines of code plots the empirical histogram:
% > [py,gridy] = empiricalHist(y);
% > figure,plot(gridy,py)

try;pr;catch;pr=1;end
[n,p] = size(y);
D = [sparse(n,1),speye(n)];
D(:,end) = [];
D = (speye(n)+D)./2;
D(n,:) = [];
kernel = exp(-pr*(-n/2:n/2).^2./n);
kernel = kernel(:)./sum(kernel);
py = zeros(n-1,p);
gridy = zeros(n-1,p);
for i=1:p
    sy = sort(y(:,i));
    [sy,ecdf] = unique(sy);
    dy = (sy(end) - sy(1))./(n-1);
    gridyi = [sy(1):dy:sy(end)]';
    ecdf = interp1(sy,ecdf./n,gridyi);
    py(:,i) = conv(diff(ecdf),kernel,'same');
    gridy(:,i) = D*gridyi;
    py(:,i) = py(:,i)./sum(py(:,i));
end
    
