function power = testPower(type,alpha,Estat,dof)
% returns the power of a test, given the expected effect size and dof
% function power = testPower(alpha,Et,dof)
% IN:
%   - type: flag for t- or F- test (can be set to 't' or 'F').
%   - alpha: statistical significance level
%   - Estat: expected effect size (z-score for a t-test, F value for an
%   F-test)
%   - dof: degrees of freedom
% OUT:
%   - power: ensuing statistical power
dof(isinf(dof)) = 1e8;
dof(dof==0) = 1e-8;
switch type
    case 't'
        tcritic = spm_invTcdf(1-alpha,dof);
        power = 1-spm_Tcdf(tcritic-Estat,dof);
    case 'F'
        if ~exist('ncfcdf')
            error('*** The statistics toolbox is required to compute F-test power. Sorry!')
        end
        Fcritic = spm_invFcdf(1-alpha,dof(1),dof(2));
        power = 1-ncfcdf(Fcritic,dof(1),dof(2),Estat);  
    otherwise
        disp('Error: this function only supports t- and F- tests!')
        power = [];
        return
end

