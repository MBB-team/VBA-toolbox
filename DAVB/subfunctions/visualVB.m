function x = visualVB(x,P,u,in)


A1 = in.A1;
A2 = in.A2;
mu3_0 = in.mu3_0;
s0 = P(1);
s1 = P(2);
s2 = P(3);
s3 = P(4);

[n1,n2] = size(A1);
n3 = size(A2,2);

mu3 = mu3_0;
Si3 = s3*eye(n3);
mu2 = A2*mu3;
Si2 = A2*Si3*A2' + s2*eye(n2);
mu1 = A1*mu2;
Si1 = A1*Si2*A1' + s1*eye(n1);

[F,PE1,PE2,PE3,PE4] = ve(mu1,Si1,mu2,Si2,mu3,Si3,u,P,in);
x = [mu1;mu2;mu3;PE1;PE2;PE3;vec(Si1);vec(Si2);vec(Si3);F;PE4];

stop = 0;
it = 1;
while ~stop
    
    Si1 = VBA_inv(eye(n1)/s0 + in.L1/s1);
    mu1 = Si1*(u./s0 + in.L1*A1*mu2/s1);
    
    Si2 = VBA_inv((A1'*in.L1*A1)/s1 + in.L2/s2);
    mu2 = Si2*(A1'*in.L1*mu1/s1 + in.L2*A2*mu3/s2);
    
    Si3 = VBA_inv((A2'*in.L2*A2)/s2 + in.L3/s3);
    mu3 = Si3*(A2'*in.L2*mu2/s2 + in.L3*mu3_0/s3);

    [F,PE1,PE2,PE3,PE4] = ve(mu1,Si1,mu2,Si2,mu3,Si3,u,P,in);
%     pause
    
    x = [x,[mu1;mu2;mu3;PE1;PE2;PE3;vec(Si1);vec(Si2);vec(Si3);F;PE4]];
    
    dF = x(end-1,end) - x(end-1,end-1);
    it = it+1;
    if abs(dF)<=in.df
        stop = 1;
    end
    
end



function [F,PE1,PE2,PE3,PE4] = ve(mu1,Si1,mu2,Si2,mu3,Si3,u,P,in)
A1 = in.A1;
A2 = in.A2;
mu3_0 = in.mu3_0;
s0 = P(1);
s1 = P(2);
s2 = P(3);
s3 = P(4);
PE1 = u-mu1;
PE2 = mu1-A1*mu2;
PE3 = mu2-A2*mu3;
PE4 = mu3-mu3_0;
F = -0.5*s0*(PE1'*PE1+trace(Si1)) ...
    -0.5*s1*(PE2'*in.L1*PE2+trace(Si2*(A1'*in.L1*A1))+trace(Si1*in.L1)) ...
    -0.5*s2*(PE3'*in.L2*PE3+trace(Si3*(A2'*in.L2*A2))+trace(Si2*in.L2)) ...
    -0.5*s3*(PE4'*in.L3*PE4+trace(Si3*in.L3)) ...
    -0.5*VBA_logDet(Si1) ...
    -0.5*VBA_logDet(Si2) ...
    -0.5*VBA_logDet(Si3);


