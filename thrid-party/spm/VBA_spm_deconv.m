function [inh] = spm_deconv(out,u,SNR)

% basic Wiener deconvolution
% function [inh] = spm_deconv(out,u,SNR)
% Let us define the above LTI convolution model:
%   out = conv(in,u) + e
% where e is a white noise process.
% This function inverts this model using Wiener deconvolution, i.e.
% estimates the input in, given the output out and the LTI convolution
% kernel u.
% IN:
%   - out: the output signal (in the time domain)
%   - u: the convolution kernel (in the time domain)
%   - SNR: the SNR level, in terms of the rqtio between the signal power
%   and the (white noise power). Attn: when SNR is not given, the
%   deconvolution operation operates a truncation of high frequencies
%   (low-pass filter at 95% of cumulative power of the frequency spectrum
%   of the convolution kernel).
% OUT:
%   - inh: the estimated input time series

out = out(:);
u = u(:);
nu = length(u);
L = length(out);
OUT = fft(out);
u = [u;zeros(L-length(u),1)];
U =  fft(u);
try % Wiener deconvolution
    G = abs(U).^2./(abs(U).^2 + SNR.^-1)./U;
    INh = OUT.*G;
catch % truncated inverse filter (when SNR is unknown)
    INh = OUT./U;
    sU = abs(U(1:L/2+1)).^2;
    scU = cumsum(sU)./sum(sU);
    i0 = find(scU>0.95);
    INh([i0;[1:L/2+1-i0(1)]'+L/2+1]) = 0;
end
inh = ifft(INh);
% inh = inh(1:L-nu);