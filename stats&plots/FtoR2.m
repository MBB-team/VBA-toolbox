function R2 = FtoR2(F,df1,df2)
R2 = 1- (1+F*df1/df2).^-1;