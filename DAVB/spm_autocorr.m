function R = spm_autocorr(y)
% computes simple autocorrelation function of signal y
% function R = spm_autocorr(y)
% Note: the autocorrelation is computed from lag 0 to lag t (length of the
% times series), and back. The symmetry of the autocorrelation function can
% be eyeballed using the following command:
% figure,plot(fftshift(R)')
% This way, the 0-lag correlation (always 1) is in the centre of the axes.
% IN:
%   - y: nXt time series matrix
% OUT:
%   - R: nXt autocorrelation matrix

[n,t] = size(y);
R = zeros(n,t*2);
% standardize y
my = mean(y,2);
sy = std(y,[],2);
sy(sy==0) = 1; % correct for identically constant time series
y = y - repmat(my,1,t);
y = diag(1./sy)*y;
% compute auto-correlation using FFT
for i=1:n
    fr = fft(y(i,:),t*2);
    tfr = fr';
    S = fr(:).*tfr(:);
    R(i,:) = ifft(S)'/t;
end
