function R2 = FtoR2(F,df1,df2)
% converts a F-test into a coefficient of determination (R2)
% function R2 = FtoR2(F,df1,df2)
% In the context of multiple regression analyses (GLMs), F-tests posses an
% interpretation in terms of percentage of variance explained (R2). This
% function derives the coefficient of determination (R2) corresponding to a
% given F-test statistic and its associated degrees-of-freedom.
% IN:
%   - F: F-test statistic
%   - df1: first F-test degree-of-freedom (usually: the smallest)
%   - df1: second F-test degree-of-freedom (usually: the biggest)
% OUT:
%   - R2: the corresponding coefficient of determination

R2 = 1- (1+F.*df1./df2).^-1;