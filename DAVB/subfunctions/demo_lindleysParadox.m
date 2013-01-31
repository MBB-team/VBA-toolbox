% this demo reproduces Lindley's paradox.

clear all

s0 = 1;
s1 = 1e-2;
m0 = -1;
m1 = 4;

n = 16;
my = -16:1e-1:16;
X = ones(n,1);
V = s1*(X*X') + s0*eye(n);
iV = VB_inv(V);

lbf = zeros(3,length(my));
M = zeros(1,length(my));

for i=1:length(my)

    dy = (my(i)-m0);
    lp0 = -(n/2)*log(2*pi) - (n/2)*log(s0) - (n/(2*s0))*dy.^2;
    dy = X*(my(i)-m1);
    lp1 = -(n/2)*log(2*pi) - (n/2)*log(s0+s1) - (1/2)*dy'*iV*dy;
    lbf(1,i) = lp1-lp0;
    
    S = 1./((1/s1)+(n/s0));
    M(i) = S*(n*my(i)/s0 + m1/s1);
    dy = my(i)-M(i);
    dm = m1-M(i);
    lp2 = -(n/2)*log(2*pi) - (n/2)*log(s0) - (n/(2*s0))*dy.^2 -(n*S/s0)/2 - VB_KL(M(i),S,m1,s1,'Normal');    
    lbf(2,i) = lp2-lp0;
    
    k1 = 1./s1;
    k0 = 1./s0;
    lbf(3,i) = -(1/2)*log((k1+n*k0)/k1) - (1/2)*(1/((1/k1)+(1/(n*k0))))*(my(i)-m1)^2 + (1/2)*n*s0*(my(i)-m0)^2;
    
    
    
end

figure,plot(my,lbf)
