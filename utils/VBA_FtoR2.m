function R2 = VBA_FtoR2 (F, df1, df2)
% // VBA toolbox //////////////////////////////////////////////////////////
%
% R2 = VBA_FtoR2 (F, df1, df2)
% converts an F-statistic into a coefficient of determination (R2)
%
% IN:
%   - F: F-test statistic
%   - df1: first F-test degree-of-freedom (usually: the smallest)
%   - df1: second F-test degree-of-freedom (usually: the biggest)
%
% OUT:
%   - R2: the corresponding coefficient of determination
%
% In the context of multiple regression analyses (GLMs), F-tests posses an
% interpretation in terms of percentage of variance explained (R2). This
% function derives the coefficient of determination (R2) corresponding to a
% given F-test statistic and its associated degrees-of-freedom.
%
% /////////////////////////////////////////////////////////////////////////

R2 = 1 - (1 + F .* df1 ./ df2) .^ (- 1);