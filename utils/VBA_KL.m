function [DKL,DJS] = VBA_KL(m1,v1,m2,v2,distrib)
% computes the KL and JS divergences between two distributions
% [DKL] = VBA_KL(m1,v1,m2,v2,disp)
% Note that the KL-divergence is not symmetric w.r.t. to the two
% probability densities, i.e. KL(p1,p2) ~=KL(p2,p1). By convention, the
% KL-divergence is defined as follows:
% KL(p1,p2) = E[log(p1)-log(p2)]
% where the expectation is taken under p1.
% IN:
%   - m1,m2: the first order moments of p1 and p2, respectively. 
%   - v1,v2: the second-orde moments of p1 and p2, respectively.
%   - distrib: the type of distribution ({'Normal'} or 'Gamma')
% OUT:
%   - DKL: the KL-divergence KL(p1,p2)
%   - DJS: the Jensen-Shannon divergence (symmetrized KL)

if isequal(m1,m2) && isequal(full(v1),full(v2))
    DKL = 0;
    DJS = 0;
    return
end

try, distrib; catch distrib='Normal'; end
DJS = [];

switch distrib

    case 'Normal'
        
        n = size(m1,1);
        iv2 = VBA_inv(v2);
        vv = v1*iv2;
        DKL = 0.5*(-VBA_logDet(vv) + trace(vv) + (m1-m2)'*iv2*(m1-m2) - n);
        
        if nargout > 1
            DJS = VBA_JensenShannon({m1,m2},{v1,v2});
        end
        
    case 'Gamma'
        
        % derive standard parameters of the Gamma distribution
        b1 = m1/v1;
        b2 = m2/v2;
        a1 = b1*m1;
        a2 = b2*m2;

        try
            psia1 = psi(a1);
        catch
            psia1 = VBA_psi(a1);
        end
        
        DKL = gammaln(a2) - gammaln(a1) - a1*(1-b2/b1) + a2*log(b1/b2) + (a1-a2)*psia1 ;
        
    otherwise
        
        DKL = [];
        
end

