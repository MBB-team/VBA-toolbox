function [F,FH0] = getF(posterior,out,Laplace)
% Free Energy lower bound on log model evidence
% function [F,FH0] = getF(posterior,out,Laplace)
% IN:
%   - posterior: the 'posterior' structure after VB inversion
%   - out: the 'out' structure after VB inversion
%   - Laplace: switch variable for the Laplace Free Energy approximation
% OUT:
%   - F: the Free Energy lower bound on the log model evidence
%   - FH0: the log evidence under the null
%------------------------------------------------------------
% Copyright (C) 2012 Jean Daunizeau / License GNU GPL v2
%------------------------------------------------------------

if nargin == 3
    out.options.Laplace = Laplace;
end

[F] = VBA_FreeEnergy(posterior,out.suffStat,out.options);

[FH0] = VBA_LMEH0(out.y);