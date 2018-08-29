function [power] = VBA_power(type, alpha, Estat, dof)
% // VBA toolbox //////////////////////////////////////////////////////////
%
% [power] = VBA_power(type, alpha, Estat, dof)
% Statistical power of a t- or F-test
%
% IN:
%   - type: type of test. Can be set to 't' or 'F'.
%   - alpha: statistical significance level
%   - Estat: expected effect size (z-score for a t-test, F value for an
%            F-test)
%   - dof: degrees of freedom
%
% OUT:
%   - power: ensuing statistical power
%
% /////////////////////////////////////////////////////////////////////////

%  ensure numerical stability
dof(isinf(dof)) = 1e8;
dof(dof==0) = 1e-8;

switch type
    case 't'
        tcritic = VBA_spm_invTcdf (1 - alpha, dof);
        power = 1 - VBA_spm_Tcdf (tcritic - Estat, dof);
    case 'F'
        Fcritic = VBA_spm_invFcdf (1 - alpha, dof(1), dof(2));
        power = 1-VBA_ncfcdf (Fcritic, dof(1), dof(2), Estat);  
    otherwise
        disp('*** VBA_power: ''type'' argument can only be ''t'' or ''F''')
        power = [];
        return
end

