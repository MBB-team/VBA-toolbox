function [g,dgdx,dgdp] = g_conv0(x,P,u,in)
% this function evaluates a linear convolution of the inputs
% function [gx,dgdx,dgdp] = g_conv0(x,P,u,in)
% This function evaluates the following multiple convolution model:
%   y(t) = cst + sum_k sum_tau w(k,tau)*u(k,t-tau)
% where u can be of any dimension (k=1,...,K).
% Note that the convolution (Volterra) kernels are constructed from the
% parameters' vector P. The size of the parameters' vector P is determined
% by the maximum lag (tau) considered.
% IN:
%   - x: [useless]
%   - P: (K*maxLag+1)x1 vector of kernel parameters.
%   - u: (K*nt)x1 vectorized input to the system
%   - in: [useless]
% OUT:
%   - g: the predicted system's output
%   - dgdx: [useless]
%   - dgdp: the gradient of the system's output w.r.t. kernel parameters
% SEE ALSO: g_convSig
if ~ isfield (in, 'K')
    in.K = 0;
end

try
    dgdp = in.dgdp;
catch
    nt = size(u,1)/in.dim.nu;
    dgdp = zeros(size(P,1),nt);
    for i=1:in.dim.nu
        ui = u((i-1)*nt+1:i*nt);
        for j=1:nt
            if j<=in.dim.n_t
                dgdp((i-1)*in.dim.n_t+j,:) = circshift(ui',[0,j-1]);
                dgdp((i-1)*in.dim.n_t+j,1:j-1) = 0;
            end
        end
    end
    dgdp(end,:) = 1;
end
g = dgdp'*P + in.K;
dgdx = [];
