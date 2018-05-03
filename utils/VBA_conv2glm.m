function [Xu] = VBA_conv2glm(u,tau)
% transforms input times series into a GLM matrix for convolution model
% function [X] = VBA_conv2glm(u,tau)
% Note that a constant term is added to the model, which means that the GLM
% convolution model contains nu*tau+1 parameters.
% IN:
%   - u: nuXnt input times series
%   - tau: the maximum lag considered in the convolution model
% OUT:
%   - Xu: (nu*tau+1)Xnt GLM matrix

try;tau;catch;tau=16;end
[nu,nt] = size(u);
np = tau*nu+1;
Xu = zeros(np,nt);
for i=1:nu
    ui = u(i,:);
    for j=1:nt
        if j<=tau
            Xu((i-1)*tau+j,:) = circshift(ui,[0,j-1]);
            Xu((i-1)*tau+j,1:j-1) = 0;
        end
    end
end
Xu(end,:) = 1;