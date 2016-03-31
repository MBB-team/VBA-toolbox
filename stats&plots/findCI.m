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
    gx = myFcdf(P,in.df(1),in.df(2));
elseif isequal(in.type,'t')
    gx = myTcdf(P,in.df);
else
    disp('findCI: error: unsupported dfistribution type!')
end

function F = myTcdf(x,v)
% Cumulative Distribution Function (CDF) of Students t distribution
% Copyright (C) 1992-2011 Wellcome Trust Centre for Neuroimaging
if nargin<2, error('Insufficient arguments'), end
ad = [ndims(x);ndims(v)];
rd = max(ad);
as = [[size(x),ones(1,rd-ad(1))];[size(v),ones(1,rd-ad(2))]];
rs = max(as);
xa = prod(as,2)>1;
if all(xa) && any(diff(as(xa,:)))
    error('non-scalar args must match in size');
end
%-Initialise result to zeros
F = zeros(rs);
%-Only defined for strictly positive v. Return NaN if undefined.
md = ( ones(size(x))  &  v>0 );
if any(~md(:))
    F(~md) = NaN;
    warning('Returning NaN for out of range arguments');
end
%-Special case: f is 0.5 when x=0 (where betainc involves log of zero)
F( md  &  x==0 ) = 0.5;
%-Special case: Standard Cauchy distribution when v=1
ml = ( md  &  v==1 ); if xa(1), mlx=ml; else mlx=1; end
F(ml) = 0.5 + atan(x(mlx))/pi;
%-Compute where defined & not special cases
Q  = find( md  &  x~=0  &  v~=1 );
if isempty(Q), return, end
if xa(1), Qx=Q; else Qx=1; end
if xa(2), Qv=Q; else Qv=1; end
%-Compute
xQxPos = x(Qx)>0;
F(Q) = xQxPos -(xQxPos*2-1).*0.5.*betainc(v(Qv)./(v(Qv)+x(Qx).^2),v(Qv)/2,1/2);


function F = myFcdf(x,v,w)
% Cumulative Distribution Function (CDF) of F (Fisher-Snedecor) distribution
% Copyright (C) 1992-2011 Wellcome Trust Centre for Neuroimaging
if nargin<2, error('Insufficient arguments'), end
%-Unpack degrees of freedom v & w from single df parameter (v)
if nargin<3
    vs = size(v);
    if prod(vs)==2
        %-DF is a 2-vector
        w = v(2); v = v(1);
    elseif vs(end)==2
        %-DF has last dimension 2 - unpack v & w
        nv = prod(vs);
        w  = reshape(v(nv/2+1:nv),vs(1:end-1));
        v  = reshape(v(1:nv/2)   ,vs(1:end-1));
    else
        error('Can''t unpack both df components from single argument')
    end
end
%-Check argument sizes
ad = [ndims(x);ndims(v);ndims(w)];
rd = max(ad);
as = [[size(x),ones(1,rd-ad(1))];[size(v),ones(1,rd-ad(2))];[size(w),ones(1,rd-ad(3))]];
rs = max(as);
xa = prod(as,2)>1;
if sum(xa)>1 && any(any(diff(as(xa,:)),1))
    error('non-scalar args must match in size'), end
%-Initialise result to zeros
F = zeros(rs);
%-Only defined for strictly positive v & w. Return NaN if undefined.
md = ( ones(size(x))  &  v>0  &  w>0 );
if any(~md(:))
    F(~md) = NaN;
    warning('Returning NaN for out of range arguments');
end
%-Non-zero where defined and x>0
Q  = find( md  &  x>0 );
if isempty(Q), return, end
if xa(1), Qx=Q; else Qx=1; end
if xa(2), Qv=Q; else Qv=1; end
if xa(3), Qw=Q; else Qw=1; end
%-Compute
F(Q) = 1 - betainc(w(Qw)./(w(Qw) + v(Qv).*x(Qx)),w(Qw)/2,v(Qv)/2);