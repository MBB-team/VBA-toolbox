% Utility
% U(x) = t1*x+t2*x^2;

close all
clear all

t2 = -0.5;
t1 = 1;

Ng =100;
x = linspace(-5,5,Ng);
U = t1.*x+t2.*(x.^2);

figure
plot(x,U);
title('Utility')


%%
X1 = x;
X2 = x;


S0 = diag([1,10]);

q = 0.25;
beta = 1;

D = zeros(Ng,Ng,2);
g = zeros(Ng,Ng);
dg = zeros(Ng,Ng,2);
for i = 1 : Ng
    
    for j = 1 : Ng
        
        x1 = X1(i);
        x2 = X2(j);
        D(i,j,1) = x1-x2;
        D(i,j,2) = x1^2-x2^2;
        dU(i,j) = [D(i,j,1),D(i,j,2)]*[t1;t2];
        g(i,j) = 1/(1+exp(-beta*(dU(i,j)) ));
        
        
        DD =[D(i,j,1);D(i,j,2)]*[D(i,j,1);D(i,j,2)]';
        
        S = q*(inv(S0) + (1-g(i,j))^2*DD )^-1 + (1-q)*( inv(S0) + g(i,j)^2*DD )^-1;
        
        crit(i,j) = trace(S);
        
    end
end


[X,Y]=meshgrid(X1,X2);

figure
meshc(X,Y,crit)

title('criteria')