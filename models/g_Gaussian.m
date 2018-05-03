function [gx] = g_Gaussian(Xt,P,ut,in)
% Gaussian convolution observation function
% function [gx] = g_Gaussian(Xt,P,ut,in)
% This function evaluates a gaussian-like observation function.
% IN:
%   - Xt: [useless]
%   - P: a 4X1 parameter vector
%   - ut: [useless]
%   - in: input structure. Contains the grid at which to evaluate the
%   gaussian bump gunction
% OUT:
%   - gx: the Gaussian bump function evaluated on the grid.

mu = P(1);
sig = exp(P(2));    % positivity constraint
A = exp(P(3));      % [id]
k = exp(P(4));
try
    grid = in.grid;
catch
    grid = -100:0.1:100;
end
grid = grid(:);
gx = A.*exp(-0.5.*(grid-mu).^2./sig);
try
    gx = gx + k.*in.input;
end
gx = convolve(gx);
gx = gx(:);


function [gx] = convolve(x)
grid = -1:0.01:1;
gauss = exp(-0.5.*grid.^2.*1e2);
gx = conv(x,gauss);

