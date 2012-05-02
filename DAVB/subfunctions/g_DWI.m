function [gx] = g_DWI(x,Phi,u,in)

% observation function for DWI B0 field data

try; in.ana; catch; in.ana = 0; end


if in.ana
    
    for i=1:9
        ind = (i-1)*6+1:i*6;
        phi = Phi(ind);
        vecA(:,i) = in.P*phi;
        miA = min(vecA(:,i));
        maA = max(vecA(:,i));
        vecA(:,i) = (vecA(:,i)-miA-1)./(maA-miA-1);
    end
    sA = reshape(vecA',3,3,[]);
    n = size(in.P,1);
    gx = zeros(n*6,1);
    for i=1:n
        indy = i:n:5*n+i;
        A = sA(:,:,i);
%         dA = sigm(diag(A));
%         ndA = reshape(sigm(vec(A),struct('G0',2,'S0',-1)),3,3);
%         A = ndA - diag(diag(ndA)) + diag(dA);
        gx(indy) = [    A(1,1).^2 + A(2,1).^2 + A(3,1).^2
            A(1,2).^2 + A(2,2).^2 + A(3,2).^2
            A(1,3).^2 + A(2,3).^2 + A(3,3).^2
            A(1,1).*A(1,2) + A(2,2).*A(2,1) + A(3,2).*A(3,1)
            A(1,1).*A(1,3) + A(2,1).*A(2,3) + A(3,3).*A(3,1)
            A(1,3).*A(1,2) + A(2,3).*A(2,2) + A(3,2).*A(3,3)    ];
    end
    
else
    n = numel(Phi)./9;
    gx = zeros(n*6,1);
    
    for i=1:n
        ind = i:n:8*n+i;
        A = reshape(Phi(ind),3,3);
        indy = i:n:5*n+i;
        gx(indy) = [    A(1,1).^2 + A(2,1).^2 + A(3,1).^2
            A(1,2).^2 + A(2,2).^2 + A(3,2).^2
            A(1,3).^2 + A(2,3).^2 + A(3,3).^2
            A(1,1).*A(1,2) + A(2,2).*A(2,1) + A(3,2).*A(3,1)
            A(1,1).*A(1,3) + A(2,1).*A(2,3) + A(3,3).*A(3,1)
            A(1,3).*A(1,2) + A(2,3).*A(2,2) + A(3,2).*A(3,3)    ];
    end
end


