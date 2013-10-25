function [k,dkdp,landmarks,dadp]=kernel_biexp(A,alpha,beta,tList)

error('not the right function')
base = (tList.^alpha).*exp(-beta*tList);
k=A.*base;

dkdp = [ base;
    k.*log(tList+eps);
    -k.*tList ];

if nargout > 2
    % time to peak: alpha/beta
    landmarks.tMax = alpha / beta ;
    landmarks.aMax = exp(-alpha).*A./((landmarks.tMax).^(-alpha));
    
    dadp=[exp(-alpha)*(alpha/beta)^alpha; 
          exp(-alpha)*log(alpha/beta)*(alpha/beta)^alpha*A;
          -(alpha*exp(-alpha)*(alpha/beta)^alpha*A)/beta];
    
end
end