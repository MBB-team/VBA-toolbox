function SNR = VBA_spm_getSNR(y,u,hrf)

% Computes fMRI SNR from controlled time windows
% function SNR = spm_getSNR(y,uu,hrf)
% IN:
%   - y: the nXt1 ROI time series
%   - u: the nuXt2 input (experimental control)
%   - hrf: the hrf (should be same sampling rate than y)
% OUT:
%   - SNR: the SNR (power ratio: signal vs noise)

ru = VBA_spm_resample(full(u),size(y,2)./size(u,2));
nu = size(u,1);
for i=1:nu
    tmp = conv(ru(i,:),hrf);
    tmp = tmp./std(tmp);
    X(:,i) = tmp(2:size(y,2)+1)';
end
iX = pinv(X'*X)*X';
ind = find(abs(mean(ru,1))>1e-1);
nreg = size(y,1);
for i=1:nreg
    beta = iX*y(i,:)';
    yc = X*beta;
    yp(:,i) = yc(ind);
    np(:,i) = y(i,ind)'-yc(ind);
end
SNR = sum(yp(:).^2)./sum(np(:).^2);