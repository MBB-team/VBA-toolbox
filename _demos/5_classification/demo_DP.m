% demo for Dirichlet Process (CRP) classification schemes

n = 1e2; % number of samples
alpha = 1; % concentration parameter of the DP

p = 6; % dimension of the data
mu0 = zeros(p,1); %1st-order moment of the prior for SS
S0 = 1e1*eye(p); %2nd-order moment of the prior for SS
% a0 = 1e2; % shape parameter
% b0 = a0*1e-1; % scale parameter


% 1- Simulate labels under CRP process
[x1] = simulate_CRP(alpha,n);
% r = 2;
% x1 = repmat([ones(1,50/r),2*ones(1,50/r)],1,r);

% 2- Draw SS of MoG components from prior
K = max(x1);
sS0 = sqrtm(S0);
m = zeros(p,K);
s = zeros(1,K);
for k=1:K
    m(:,k) = mu0 + sS0*randn(p,1);
    s(k) = 1e-2;
end


% 2- Simulate data sample from CRP labels
y = zeros(p,n);
for i=1:n
    y(:,i) = m(:,x1(i)) + sqrt(s(x1(i)))*randn(p,1);
end

[post,out] = VB_CRP(y,alpha,struct('verbose',0));

hf = figure;
subplot(2,2,1),plot(1:n,x1),title('simulated CRP labels')
subplot(2,2,2),plot(1:n,out.xih),title('estimated CRP labels')
subplot(2,2,3),imagesc(m),title('simulated component means'),colorbar
subplot(2,2,4),imagesc(post.mu),title('estimated component means'),colorbar


