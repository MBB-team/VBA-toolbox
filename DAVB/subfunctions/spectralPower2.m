function [gy,P,gridw,g,R] = spectralPower2(P,gridw,g,C)


nw = numel(gridw); % number of frequency bins
nk = 32; % number of harmonics in the decomposition
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

g = g.*v*c.*a/2;

try
    C;
catch
    C = [1 0 0];
end

J0 = zeros(nk,1);

for w=1:nw
    
    W = gridw(w);
    
    for k=1:nk
        
        K = k.*pi./P.l;

        J0(k) = -v.*2*c.*(1-K./c.^2);
         
%         J = [   0               1           0
%                 -P.Ke.^2        -2*P.Ke     P.me.*P.Ke
%                 g               0           J0(k)          ];
%         
%         R(:,k,w) = inv(1i.*W*eye(3) - J)*ones(3,1);

        
        a11 = 1i.*W;
        a22 = 1i.*W + 2.*P.Ke;
        a33 = 1i.*W - J0(k);
        
        detF = a33.*(a11*a22+ P.Ke.^2) - P.me.*P.Ke.*g;
          
        R(1,k,w) = ( a33.*(a22+1) + P.me.*P.Ke )./detF;
        
        R(2,k,w) = ( P.me.*P.Ke.*(g+a11) + (a11-P.Ke.^2).*a33 )./detF;
        
        R(3,k,w) = ( a11.*a22 - P.Ke.^2 + (a22+1).*g )./detF;
        
        
        [L] = leadField(K,P);
        
        gy(w) = gy(w) + C*R(:,k,w)*L;
%         L.*T(k,w).*conj(T(k,w)).*conj(L);
        
    end
%     
%     J0
%     pause
    
end




function [L] = leadField(K,P)
L = P.phi(1).*exp(-P.phi(2).*pi.^2.*K.^2);

function [g] = dsdv(P)
s0 = 1./(1+exp(P.sig.r.*P.sig.eta));
g = P.sig.r.*s0.*(1-s0);

