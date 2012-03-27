function [fx,dF_dX,dF_dTheta,d2F_dXdTheta] = f_SHC(Xt,Theta,ut,inF)
% stable heteroclinic channels evolution function
%------------------------------------------------------------
% Copyright (C) 2012 Jean Daunizeau / License GNU GPL v2
%------------------------------------------------------------
deltat = inF.deltat;
x = Xt;
try
    K = inF.K;
catch
    K = [1/2;1/32];
end
try
    lambda = inF.lambda;
catch
    lambda = 0.3;
end
try
    G0 = inF.G0;
catch
    G0 = 50;
end
try
    beta = inF.beta;
catch
    beta = 0.5;
end
try
    ind1 = inF.ind1;
    R1   = inF.R1;
catch
    ind1 = 1:4;
    R1{1} = [1  5  5  .5
            .5  1  5  5
             5  .5 1  5
             5  5  .5 1];
    R1{2} = [1  5  5  .5
             5  1  .5 5
             .5 5  1  5
             5  .5 5  1];
    R1{3} = [1 .5  5  5
             5  1  .5 5
             5  5  1  .5
             .5 5  1  5];
end
try
    ind2 = inF.ind2;
    R2   = inF.R{2};
catch
    ind2 = 5:7;
    R2   = [1   5 .5
            .5  1  5
            5  .5  1 ];
end


% Separate states
x1 = x(ind1);
x2 = x(ind2);

% High level
[SX2,dsdx2] = sigm(x2,G0,beta);
ff{2} = K(2).*(-lambda*x2 - R2*SX2);

% Low level
[SX1,dsdx1] = sigm(x1,G0,beta);
R1eff = R1{1}.*SX2(1) + R1{2}.*SX2(2) + R1{3}.*SX2(3);
ff{1} = K(1).*(-lambda*x1 - R1eff*SX1);

f = [ff{1};ff{2}];

df1 = zeros(size(f,1),length(ind1));
df2 = zeros(size(f,1),length(ind2));

df2(ind2,:) = K(2).*lambda.*(-eye(length(ind2))) -K(2).*R2*diag(dsdx2);
df1(ind1,:) = K(1).*lambda.*(-eye(length(ind1))) -K(1).*R1eff*diag(dsdx1);
% df1(ind2(1),ind1) = [-K(1).*dsdx2(1).*R1{1}*SX1]';
% df1(ind2(2),ind1) = [-K(1).*dsdx2(2).*R1{2}*SX1]';
% df1(ind2(3),ind1) = [-K(1).*dsdx2(3).*R1{3}*SX1]';

df_dx = [df1,df2]';
df_dx(ind2(1),ind1) = [-K(1).*dsdx2(1).*R1{1}*SX1]';
df_dx(ind2(2),ind1) = [-K(1).*dsdx2(2).*R1{2}*SX1]';
df_dx(ind2(3),ind1) = [-K(1).*dsdx2(3).*R1{3}*SX1]';

fx = Xt + deltat.*f;
dF_dX = eye(length(Xt)) + deltat.*df_dx;
dF_dTheta = [];
d2F_dXdTheta = [];

% dfdx2 = numericDiff(@get_fx,1,x,Theta,deltat);
% 
% hf=figure,plot(dfdx2(:),dF_dX(:),'.')
% pause
% close(hf)


function [Sx,dsdx] = sigm(x,G0,beta)
bx = beta*x;
indn = find(bx<=-1e2);
indp = find(bx>=1e2);
Sx = G0./(1+exp(-bx));
dsdx = G0./(1+exp(-bx)).^2.*beta.*exp(-bx);
Sx(indn) = 0;
Sx(indp) = G0;
dsdx(indn) = 0;
dsdx(indp) = 0;


% 
% function [fx] = get_fx(x,theta,deltat)
% try
%     K = theta.K;
% catch
%     K = [1/2;1/32];
% end
% try
%     lambda = theta.lambda;
% catch
%     lambda = 0.3;
% end
% try
%     G0 = theta.G0;
% catch
%     G0 = 50;
% end
% try
%     beta = theta.beta;
% catch
%     beta = 0.5;
% end
% try
%     ind1 = theta.ind1;
%     R1   = theta.R1;
% catch
%     ind1 = 1:4;
%     R1{1} = [1  5  5  .5
%             .5  1  5  5
%              5  .5 1  5
%              5  5  .5 1];
%     R1{2} = [1  5  5  .5
%              5  1  .5 5
%              .5 5  1  5
%              5  .5 5  1];
%     R1{3} = [1 .5  5  5
%              5  1  .5 5
%              5  5  1  .5
%              .5 5  1  5];
% end
% try
%     ind2 = theta.ind2;
%     R2   = theta.R{2};
% catch
%     ind2 = 5:7;
%     R2   = [1   5 .5
%             .5  1  5
%             5  .5  1 ];
% end
% 
% 
% % Separate states
% x1 = x(ind1);
% x2 = x(ind2);
% 
% % High level
% [SX2,dsdx2] = sigm(x2,G0,beta);
% ff{2} = K(2).*(-lambda*x2 - R2*SX2);
% 
% % Low level
% [SX1,dsdx1] = sigm(x1,G0,beta);
% R1eff = R1{1}.*SX2(1) + R1{2}.*SX2(2) + R1{3}.*SX2(3);
% ff{1} = K(1).*(-lambda*x1 - R1eff*SX1);
% 
% f = [ff{1};ff{2}];
% 
% fx = x + deltat.*f;
