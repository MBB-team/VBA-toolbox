function [gy,P,gridw,g,T] = spectralPower(P,gridw,g)
%------------------------------------------------------------
% Copyright (C) 2012 Jean Daunizeau / License GNU GPL v2
%------------------------------------------------------------

nw = numel(gridw); % number of frequency bins
nk = 32; % number of harmonics in the decomposition
gy = zeros(1,nw);
R =  zeros(nk,nw);
T =  zeros(nk,nw);
D = zeros(3,3,nk,nw);

try % synaptic gain at steady-state membrane depolarization
    g;
catch
    [g] = dsdv(P);
end

for w=1:nw
    W = pi.*gridw(w);
    
    for k=1:nk
        K = k.*pi./P.l;
        
        D(:,:,k,w) = Dij(K,W,P);
        
        R(k,w) = (P.Ki+1i*W).^2.*( ...
            P.Ke.^4 + 4*1i*P.Ke.^3.*W - 4*1i*P.Ke.*W.^3 + W.^4 - P.Ke.^2.* ...
            (D(1,3,k,w).*D(3,1,k,w).*g.^2.*P.me.^2 + 6.*W.^2) ) ...
            -D(2,3,k,w).*D(3,2,k,w).*g.^2.*P.Ke.*P.Ki.*P.me.*P.mi.*(P.Ke+1i*W).^2;
        
        T(k,w) = D(3,1,k,w).*g.*P.Ke.^2.*P.me.^2.*(P.Ki+1i).^2./R(k,w);
        
        [L] = leadField(K,P);
        
        gy(w) = gy(w) + (pi./P.l).*(real(T(k,w)).^2+imag(T(k,w)).^2).*L.^2;
%         L.*T(k,w).*conj(T(k,w)).*conj(L);
        
    end
end




function [D] = Dij(K,W,P)
D = zeros(3,3);
ind = {[1,3],[3,1],[2,3],[3,2]};
for i=1:length(ind)
    aij = P.a(ind{i}(1),ind{i}(2));
    cij = P.c(ind{i}(1),ind{i}(2));
    vij = 1./P.v(ind{i}(1),ind{i}(2));
    D(ind{i}(1),ind{i}(2)) = aij.*(cij+1i*vij*W)...
        ./(cij.^2 - vij.^2.*W.^2 + 2*1i*vij.*cij.*W + K.^2);
end


function [L] = leadField(K,P)
L = P.phi(1).*exp(-P.phi(2).*pi.^2.*K.^2);

function [g] = dsdv(P)
g = P.sig.r.*exp(P.sig.r.*P.sig.eta)./(1+exp(P.sig.r.*P.sig.eta)).^2;

