function [Sx,dsdx,dsdp] = sigm(x,in,Phi)
% evaluates the sigmoid function (weird parameterization)
% function [Sx,dsdx,dsdp] = sigm(x,in,Phi)
% IN:
%   - x: the value at which the sigmoid function is evaluated. Can be a
%   vector, in which case the output is also a (row) vector.
%   - in: a structure that contains some default gains:
%       .G0: is the factor by which the sigmoid function is scaled ({1})
%       .S0: is the intercept ({0})
%       .beta: is the default slope of the sigmoid ({1})
%       .INV: if 1, the function evaluates the inverse sigmoid function
%       (this means sigm(sigm(x,struct('INV',1))) = x)
%       .mat: if 1, do not vectorize x
%   - Phi: A vector, can be left unspecified, or contain one or two
%   entries. The first entry rescale the slope in.beta ({0}); the second is
%   the inflexion point ({0})
% OUT:
%   - Sx: the sigmoid function evaluated at x
%   - dsdx: the derivative of the sigmoid function wrt x
%   - dsdp: the derivative of the sigmoid function wrt Phi.

% get default parameterization
if ~exist('in','var')
    in = [];
end
in = VBA_check_struct(in,'G0',1,'S0',0,'beta',1,'INV',0,'deriv',0,'mat',0);
if ~in.mat
    % transform x into a row vector
    x = x(:)';
end

if exist('Phi','var') && ~isempty(Phi)
    b = in.beta.*exp(Phi(1));
    if size(Phi,1) > 1
        th = Phi(2);
    else
        th = 0;
    end
else
    b = in.beta;
    Phi = [];
    th = 0;
end

if in.INV
    % evaluate inverse sigmoid
    ixg = in.G0.*(x-in.S0).^-1 - 1;
    Sx = th - b.^-1.*log(ixg);
    dsdx = [];
    dsdp = [];
    return
else
    % evaluate sigmoid (and check numerical instabilities)
    bx = b*(x-th);
    Sx = in.S0 + in.G0./(1+exp(-bx));
end

% Sx(Sx<1e-3)=1e-3;
% Sx(Sx>1-1e-3)=1-1e-3;

if nargout < 2 ; return; end

% evaluate derivative wrt x
dsdx = b*(Sx-in.S0).*(1-(Sx-in.S0)./in.G0);

if nargout < 3 ; return; end


% evaluate derivatives wrt parameters
if size(Phi,1) >= 1
    dsdp = zeros(size(Phi,1),length(x));
    dsdp(1,:) = (x-th).*in.beta.*dsdx;
    if size(Phi,1) == 2
        dsdp(2,:) = -dsdx;
    end
    dsdp(isnan(dsdp)) = 0;
else
    dsdp = [];
end
