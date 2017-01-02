function Z = RFT_Gtf(t,nu)
% derives the Gaussianized Z RF from t-field with nu d.o.f.
% function Z = RFT_Gtf(t,nu)
% IN:
%   - t: t-score
%   - nu: d.o.f. of the t-score
% OUT:
%   - Z: Gaussianized t-score
Z =  icdf('norm',tcdf(t,nu));