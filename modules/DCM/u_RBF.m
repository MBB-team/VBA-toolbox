function [u,dudx,dudp] = u_RBF(x,P,t,in)
% input function for free-form (RBF) deterministic DCM
% function [u,dudx,dudp] = u_RBF(x,P,t,in)
% IN:
%   - x: [useless]
%   - P: nW*nuX1 vector of RBF projection parameters, where nW is the
%   number of RBFs and nu is the dimensionality of the output
%   - t: time at which the RBF set is evaluated
%   - in: user-defined structure containing the set of RBF parameters
%   (centres and fwhm for gaussian RBFs)
% OUT:
%   - u: the nuXnt matrix of outputs, where nt is the length of the
%   argument t.
%   - dudx: [useless]
%   - dudp: gradient wrt the parameters

% construct Fourier basis at time t
X = zeros(length(in.centres),length(t));
for i = 1:length(t)
    X(:,i) = exp(-0.5*(in.centres(:)-t(i)).^2./in.sig)./sqrt(2*pi*in.sig);
end

% project basis onto input space
P = reshape(P,length(in.centres),[]);
u = (X'*P)';

% get gradients
dudx = [];
dudp = kron(eye(size(P,2)),X);