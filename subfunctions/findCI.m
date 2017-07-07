function [ci] = findCI(alpha,type,dof)
% finds confidence intervals for given t or F distributions
% function [ic] = findCI(alpha,dof,type)
% This function uses VBA to find the values x1 and x2, such that:
% P(x>x1) = alpha and P(x>x1) = 1-alpha, where P is the cumulative density
% function of either t- or F- distributions and alpha is the confidence
% level.
% IN:
%   - alpha: confidence level (0<alpha<1)
%   - type: distribution type ('t' or 'F')
%   - dof: corresponding degrees of freedom. NB: for F-distributions, dof
%   is a 2x1 vector!
% OUT:
%   - ci: confidence interval

options.priors.SigmaPhi = 1e4;
options.inG.type = type;
options.inG.df = dof;
options.verbose = 0;
options.DisplayWin = 0;
dim.n_phi = 1;
dim.n = 0;
dim.n_theta = 0;
[post] = VBA_NLStateSpaceModel(alpha,[],[],@pval,dim,options);
ci(1) = post.muPhi;
[post] = VBA_NLStateSpaceModel(1-alpha,[],[],@pval,dim,options);
ci(2) = post.muPhi;

function gx = pval(x,P,u,in)
if isequal(in.type,'F')
    gx = spm_Fcdf(P,in.df(1),in.df(2));
elseif isequal(in.type,'t')
    gx = spm_Tcdf(P,in.df);
else
    disp('findCI: error: unsupported dfistribution type!')
end