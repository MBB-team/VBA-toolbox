function [gy,P,gridw,g,R] = spectralPower3(P,gridw,g,Cout)


nw = numel(gridw); % number of frequency bins
nk = 16; % number of harmonics in the decomposition
gy = zeros(1,nw);
R =  zeros(3,nk,nw);

try % synaptic gain at steady-state membrane depolarization
    g;
catch
    [g] = dsdv(P);
end
v = P.v(P.i1,P.i2);
c = P.c(P.i1,P.i2);
a = P.a(P.i1,P.i2);

cz = 0.5*g.*v.*a/c;


try
    Cout;
catch
    Cout = [1 0 0];
end

Cin = [0;1;0];


J0 = zeros(nk,1);

for w=1:nw
    
    W = gridw(w)*2*pi;
    
    for k=1:nk
        
        K = k*pi/P.l;
        fk = -0.5*v*(1-K*c^2)/c;

        J = [   0           0           1
                cz          fk          0
                -P.Ke^2     P.me*P.Ke   -2*P.Ke  ];
            
        R = VB_inv(1i.*W*eye(3) - J)*Cin;

        L=1;%[L] = leadField(K,P);
        
        gy(w) = gy(w) + Cout*R*L;
%         L.*T(k,w).*conj(T(k,w)).*conj(L);
        
    end

    
end




function [L] = leadField(K,P)
L = P.phi(1).*exp(-pi.^2./(P.phi(2).*K.^2));

function [g] = dsdv(P)
s0 = 1./(1+exp(P.sig.r.*P.sig.eta));
g = P.sig.r.*s0.*(1-s0);

