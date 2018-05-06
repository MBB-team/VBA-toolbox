function [sg,dsdx,dgdp] = g_convSig(x,P,u,in)
% this function evaluates a logistic convolution of the inputs
% function [gx,dgdx,dgdp] = g_convSig(x,P,u,in)
% This function evaluates the following multiple convolution model:
%   p(y(t)=1) = sigm(cst + sum_k sum_tau w(k,tau)*u(k,t-tau))
% where u can be of any dimension (k=1,...,K), and sigm is the standard
% sigmoid mapping. This logistic convolution model is useful for binary
% observations.
% Note that the convolution (Volterra) kernels are constructed from the
% parameters' vector P. The size of the parameters' vector P is determined
% by the maximum lag (tau) considered.
% IN:
%   - x: [useless]
%   - P: (K*maxLag+1)x1 vector of kernel parameters.
%   - u: (K*nt)x1 vectorized input to the system
%   - in: [useless]
% OUT:
%   - sg: the predicted system's output (in terms of p(y(t)=1))
%   - dsdx: [useless]
%   - dsdp: the gradient of the system's output w.r.t. kernel parameters
% SEE ALSO: g_conv0

[g,dsdx,dgdp] = g_conv0(x,P,u,in);
sg = VBA_sigmoid(g);
dgdp = dgdp*diag(sg.*(1-sg));