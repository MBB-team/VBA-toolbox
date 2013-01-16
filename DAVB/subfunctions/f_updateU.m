function fx = f_updateU(x,P,u,in)

% choices
ac = u(in.ic).*u(in.iu(1)) + (1-u(in.ic)).*u(in.iu(2));
au = u(in.ic).*u(in.iu(2)) + (1-u(in.ic)).*u(in.iu(1));
% pause
% u(:,t) = c*ind(:) + (1-c)*flipud(ind(:));
% ac = u(in.ic); % chosen item
% au = u(in.iu); % unchosen item
% previous belief about utility
mu0 = x(1:in.n); % expectation (E[x] = mu)
S0 = reshape(x(in.n+1:in.n+in.n^2),in.n,in.n); % variance (V[x] = mu)
% intermediary variables
dV = x(ac) - x(au); % value of chosen item (relative to unchosen item)
W = eye(in.n); % utility basis function set
ddVdx = W(:,ac)-W(:,au); % gradient of dV
b = exp(P(in.temp)); % behavioural temperature
g = sig(dV/b); % probability of picking the chosen item
% VB update rule
S = pinv(pinv(S0) + (1-g)*g*(ddVdx*ddVdx'));
mu = mu0 + S*ddVdx*(1-g);
% wrap up
fx = zeros(size(x));
fx(1:in.n) = mu;
fx(in.n+1:in.n+in.n^2) = vec(S);


function s= sig(x)
s = 1./(1+exp(-x));
s(s<1e-3) = 1e-3;
s(s>1-1e-3) = 1-1e-3;