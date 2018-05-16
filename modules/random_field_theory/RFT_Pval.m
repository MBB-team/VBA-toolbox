function [Pw] = RFT_Pval(u,k,c,fwhm,L,type,dof)
% computes Pval according to Friston et al. (1996), Eq. 1-4.
% [Pw] = RFT_Pval(u,k,c,fwhm,L)
% This function ca be called to derive the p-value at different levels of
% inference, i.e.:
% - P_peak = RFT_Pval(u,0,1,fwhm,L)
% - P_cluster = RFT_Pval(U,k,1,fwhm,L), where U was the set-inducing
% threshold, and k is the observed spatial extent of the cluster.
% - P_set = RFT_Pval(U,K,c,fwhm,L), where U and K were the set-inducing
% thresholds, and c was the number of observed upcrossing clusters
% IN:
%   - u: RF value
%   - k: cluster's spatial extent
%   - c: number of upcrossing clusters
%   - fwhm: estimated FWHM
%   - L: search volume
%   - type: type of RF. Can be set to 'norm' (normal, default), 't'
%   (Student) or 'F' (Fisher).
%   - dof: degrees of freedom (only relevant for 't' or 'F' fields).
% OUT:
%   - Pw: the ensuing p-value

EC = RFT_expectedTopo(u,L,fwhm,1,type,dof);
switch type
    case 'norm'
        P0 = 1-VBA_spm_Ncdf(u,0,1);
    case 't'
        P0 = 1-VBA_spm_Tcdf(u,dof);
    case 'F'
        P0 = 1-VBA_spm_Fcdf(u,dof(1),dof(2));
end
beta = (gamma(3/2).*EC./(L.*P0)).^2;
Pnk = exp(-beta.*k.^2);
if k == 0, Pnk = 1; end; % solves the numerical issue when P0=0    
Pw = 1;
for i=0:c-1
    Pw = Pw - myPoissonPMF(i,EC.*Pnk);
end

function p = myPoissonPMF(x,Ex)
p = (Ex.^x).*exp(-Ex)./factorial(x);

% function p = myPoissonCDF(x,Ex)
% p = 1 - gammainc(Ex,x+1);


