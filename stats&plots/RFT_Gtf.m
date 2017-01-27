function Z = RFT_Gtf(x,df,type)
% derives the Gaussianized Z RF from t-field with nu d.o.f.
% function Z = RFT_Gtf(x,dof,type)
% IN:
%   - x: t-score or F-score
%   - nu: d.o.f.
%   - type: 't' or 'F'
% OUT:
%   - Z: Gaussianized z-score
try, type; catch, type = 't'; end
switch type
    case 't'
        Z =  icdf('norm',tcdf(x,df),0,1);
    case 'F'
        Z =  icdf('norm',fcdf(x,df(1),df(2)),0,1);
    otherwise
        disp('Error: RFT_Gtf: type is either ''t'' or ''F''!')
        Z = [];
end