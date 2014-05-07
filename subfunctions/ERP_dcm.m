function [fx,dF_dX,dF_dTheta] = ERP_dcm(x,Theta,u,in)
% neural mass evolution function (DCM for ERPs)
% function [fx,dF_dX,dF_dTheta,d2F_dXdTheta] =
% ERP_dcm(Xt,Theta,ut,in)
%
% This function ouputs the state-transition associated with the
% 5-populations neural mass microcolumn:
%   - 2 inhibitory populations in the supragranular layers
%   - 1 excitatory population in the granular layers
%   - 2 excitatory populations (pyramidal cells) in the infragranular
%   layers.

delta_t = in.delta_t;

n = size(x,1);

He = Theta(1);
Hi = Theta(2);
Ke = Theta(3);
Ki = Theta(4);
Ka = Theta(5);
g1 = Theta(6);
g2 = Theta(7);
g3 = Theta(8);
g4 = Theta(9);
g5 = Theta(10);
r1 = Theta(11);
r2 = Theta(12);
a = Theta(13);

[sx,g] = sig(x,r1,r2);


xDot = [    x(4)
            x(5)
            x(6)
            Ke.*He.*(g1.*sx(9)+u)-2.*Ke.*x(4)-Ke.^2.*x(1)
            Ke.*He.*g2.*sx(1)-2.*Ke.*x(5)-Ke.^2.*x(2)
            Ki.*Hi.*g4.*sx(12)-2.*Ki.*x(6)-Ki.^2.*x(3)
            x(8)
            Ke.*He.*g3.*sx(9)-2.*Ke.*x(8)-Ke.^2.*x(7)
            x(5)-x(6)
            x(11)
            Ki.*Hi.*g5.*sx(12)-2.*Ki.*x(11)-Ki.^2.*x(10)
            x(8)-x(11)];
        




A = [0                  0       0       1       0       0       0       0       0                   0       0       0
     0                  0       0       0       1       0       0       0       0                   0       0       0
     0                  0       0       0       0       1       0       0       0                   0       0       0
     -Ke.^2             0       0       -2.*Ke  0       0       0       0       Ke.*He.*g1.*g(9)    0       0       0
     Ke.*He.*g2.*g(1)   -Ke.^2  0       0       -2.*Ke  0       0       0       0                   0       0       0
     0                  0       -Ki.^2  0       0       -2.*Ki  0       0       0                   0       0       Ki.*Hi.*g4.*g(12)
     0                  0       0       0       0       0       0       1       0                   0       0       0
     0                  0       0       0       0       0       -Ke.^2  -2.*Ke  Ke.*He.*g3.*g(9)    0       0       0
     0                  0       0       0       1       -1      0       0       0                   0       0       0
     0                  0       0       0       0       0       0       0       0                   0       1       0
     0                  0       0       0       0       0       0       0       0                   -Ki.^2  -2.*Ki  Ki.*Hi.*g5.*g(12)
     0                  0       0       0       0       0       0       1       0                   0       -1       0                   ];


 
     
B = [  [0,   0,   0,   Ke.*(g1.*sx(9)+u),                          Ke.*g2.*sx(1),                           0,                                      0,   Ke.*g3.*sx(9),                        0,   0,   0,                                         0 ];
       [0,   0,   0,   0,                                          0,                                       Ki.*g4.*sx(12),                         0,   0,                                    0,   0,   Ki.*g5.*sx(12),                            0 ];
       [0,   0,   0,   He.*(g1.*sx(9)+u)-2.*x(4)-2.*Ke.*x(1),      He.*g2.*sx(1)-2.*x(5)-2.*Ke.*x(2),       0,                                      0,   He.*g3.*sx(9)-2.*x(8)-2.*Ke.*x(7),    0,   0,   0,                                         0 ];
       [0,   0,   0,   0,                                          0,                                       Hi.*g4.*sx(12)-2.*x(6)-2.*Ki.*x(3),     0,   0,                                    0,   0,   Hi.*g5.*sx(12)-2.*x(11)-2.*Ki.*x(10),      0 ];
       [0,   0,   0,   0,                                          0,                                       0,                                      0,   0,                                    0,   0,   0,                                         0 ];
       [0,   0,   0,   Ke.*He.*sx(9),                              0,                                       0,                                      0,   0,                                    0,   0,   0,                                         0 ];
       [0,   0,   0,   0,                                          Ke.*He.*sx(1),                           0,                                      0,   0,                                    0,   0,   0,                                         0 ];
       [0,   0,   0,   0,                                          0,                                       0,                                      0,   Ke.*He.*sx(9),                        0,   0,   0,                                         0 ];
       [0,   0,   0,   0,                                          0,                                       Ki.*Hi.*sx(12),                         0,   0,                                    0,   0,   0,                                         0 ];
       [0,   0,   0,   0,                                          0,                                       0,                                      0,   0,                                    0,   0,   Ki.*Hi.*sx(12),                            0 ];
       [0,   0,   0,   0,                                          0,                                       0,                                      0,   0,                                    0,   0,   0,                                         0 ];
       [0,   0,   0,   0,                                          0,                                       0,                                      0,   0,                                    0,   0,   0,                                         0 ];
       [0,   0,   0,   0,                                          0,                                       0,                                      0,   0,                                    0,   0,   0,                                         0 ]    ];
     

     
     
% state transition (Euler discretization)
fx = x + delta_t*xDot;

% Jacobian
dF_dX = (speye(n) + delta_t*A)';


% gradient wrt the evolution parameters
dF_dTheta = delta_t*B;



