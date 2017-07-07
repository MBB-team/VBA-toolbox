function [h,p,KSstatistic]=VBA_kstest(x)
% alternative to kstest fromt the statistics toolbox

%% compute critical statistic
n = numel(x);

sampleCDF = linspace(0,1,n+1)';
x = sort(x);

nullCDF  =  spm_Ncdf(x,0,1);

delta1    =  sampleCDF(1:end-1) - nullCDF;   % Vertical difference at jumps approaching from the LEFT.
delta2    =  sampleCDF(2:end)   - nullCDF;   % Vertical difference at jumps approaching from the RIGHT.
deltaCDF  =  abs([delta1 ; delta2]);

KSstatistic   =  max(deltaCDF);

%% compute p-value
s = n*KSstatistic^2;

% For d values that are in the far tail of the distribution (i.e.
% p-values > .999), the following lines will speed up the computation
% significantly, and provide accuracy up to 7 digits.
if (s > 7.24) ||((s > 3.76) && (n > 99))
    p = 2*exp(-(2.000071+.331/sqrt(n)+1.409/n)*s);
else
    % Express d as d = (k-h)/n, where k is a +ve integer and 0 < h < 1.
    k = ceil(KSstatistic*n);
    h = k - KSstatistic*n;
    m = 2*k-1;
    
    % Create the H matrix, which describes the CDF, as described in Marsaglia,
    % et al.
    if m > 1
        c = 1./gamma((1:m)' + 1);
        
        r = zeros(1,m);
        r(1) = 1;
        r(2) = 1;
        
        T = toeplitz(c,r);
        
        T(:,1) = T(:,1) - (h.^(1:m)')./gamma((1:m)' + 1);
        
        T(m,:) = fliplr(T(:,1)');
        T(m,1) = (1 - 2*h^m + max(0,2*h-1)^m)/gamma(m+1);
    else
        T = (1 - 2*h^m + max(0,2*h-1)^m)/gamma(m+1);
    end
    
    % Scaling before raising the matrix to a power
    if ~isscalar(T)
        lmax = max(eig(T));
        T = (T./lmax)^n;
    else
        lmax = 1;
    end
    
    % Pr(Dn < d) = n!/n * tkk ,  where tkk is the kth element of Tn = T^n.
    % p-value = Pr(Dn > d) = 1-Pr(Dn < d)
    p = (1 - exp(gammaln(n+1) + n*log(lmax) - n*log(n)) * T(k,k));
end

% test hypothesis
h  =  (p < .05);

end
