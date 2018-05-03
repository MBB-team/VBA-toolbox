function [C] = Contrast_MEbins(n)
% F-contrast for a main effect of a VI with n different values
C = zeros(n,n*(n-1)/2);
ii = 1;
for i=1:n-1
    for j=i+1:n
       C(i,ii) = 1;
       C(j,ii) = -1;
       ii = ii + 1;
    end
end

