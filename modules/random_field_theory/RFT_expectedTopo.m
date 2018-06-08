function [HC,EC] = RFT_expectedTopo(u,L,fwhm,c,type,dof)
% derives E[EC] and E[HC] for a Gaussian RF defined on a regular lattice
% [out] = RFT_expectedTopo(u,L,fwhm,c,verbose)
% IN:
%   - u: set-inducing threshold
%   - L: length of the lattice
%   - fwhm: smoothness of the random field
%   - c: number of unbroken field segments (default=1)
%   - type: type of RF. Can be set to 'norm' (normal, default), 't'
%   (Student) or 'F' (Fisher).
% OUT:
%   - HC: expected hadwiger characteristic
%   - EC: expected Euler characteristics

try;c;catch;c=1;end
try,type;catch;type='norm';end
try,dof;catch;dof=NaN;end
u = VBA_vec(u);
switch type
    case 'norm'
        P = 1-VBA_spm_Ncdf(u,0,1);
        EC = (L./fwhm).*(sqrt(4*log(2))./(2*pi)).*exp(-u.^2./2);
        HC = c*P + EC;
    case 't'
        P = 1 - VBA_spm_Tcdf(u,dof);
        EC = (L./fwhm).*(sqrt(4*log(2))/(2*pi))*(1+u.^2/dof).^((1-dof)/2);
        HC = c*P + EC;
    case 'F'
        b = gammaln(dof(2)/2) + gammaln(dof(1)/2);
        P = 1 - VBA_spm_Fcdf(u,dof(1),dof(2));
        EC = (L./fwhm).*(sqrt(4*log(2))/(2*pi)).*exp(gammaln((sum(dof)-1)/2)-b)*2^(1/2)...
            *(dof(1)*u/dof(2)).^(1/2*(dof(1)-1)).*(1+dof(1)*u/dof(2)).^(-1/2*(sum(dof)-2));
        HC = c*P + EC;
    otherwise
        disp('Error: RFT_expectedTopo: unknown type. Please use ''norm'', ''t'' or ''F''.')
        HC = [];
        EC = [];
end

