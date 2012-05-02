% demo lag


lag = 5;
n = 2;
p = 5;


g = randn(p,1);
y = randn(p,1);

f = randn(n,1);

G = randn(p,n);
F = randn(n,n);

mutm1 = randn(n*lag,1);
Stm1 = randn(n*lag,n*lag);
Stm1 = Stm1'*Stm1;

m = randn(n*lag,1);

Q = eye(p);
R = eye(n);


X = randn(n*lag,1);

costFcn = @lagged_p;
init = 0*randn(size(m));
options.args = {g,y,Q,G,R,F,f,mutm1,Stm1,m,lag};
options.minimize = 0;

nit = 1;
xx = zeros(n*lag,nit);
for i=1:nit
    init = 0.*randn(n*lag,1);
    [mu,curv,out] = optimCost(costFcn,init,options);
    xx(:,i) = mu;
end





% numericDiff(@numericDiff,3,costFcn,1,X,options.args{:})
% 
% 
[Q,St,mut] = lagged_p(mu,g,y,Q,G,R,F,f,mutm1,Stm1,m,lag);
dcdx = numericDiff(costFcn,1,mu,options.args{:})
dcdx2 = numericDiff(costFcn,1,mut,options.args{:})
mut - mu
-curv - St

