function [x1] = simulate_CRP(alpha,n)
% simulate labels from CRP process with concentration alpha
% IN:
%   - alpha: concentration parameter of the CRP
%   - n: size of the time series
% OUT:
%   - x1: CRP-label process (vector of size 1Xn)

% initial condition
x1 = zeros(1,n);
x1(1) = 1;

% loop over time samples
for i=2:n
   Ki = max(x1(1:i-1)); % number of classes so far
   nk = zeros(Ki,1);
   for k=1:Ki
       nk(k) = sum(x1(1:i-1)==k);
   end
   p = [nk;alpha]./(alpha+i-1);
   x1(i) = VBA_random ('Categorical', p);
   
end