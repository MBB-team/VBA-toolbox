function [k,dkdp,landmarks,dadp]=kernel_biexp(A,alpha,beta,tList)

% error('*** bad gradients');

tList = tList+2*eps;
k=A.*(tList.^alpha).*exp(-beta*tList);

dkdp = [ (tList.^alpha).*exp(-beta*tList);
    	 A*(tList.^alpha).*exp(-beta*tList).*log(tList)
        -A*(tList.^(alpha+1)).*exp(-beta*tList) ]; 

if nargout > 2
    % time to peak: alpha/beta
    landmarks.tMax = alpha / beta ;
    landmarks.aMax = exp(-alpha).*A./((landmarks.tMax).^(-alpha));
    
    dadp=[exp(-alpha)*(alpha/beta)^alpha; 
          exp(-alpha)*log(alpha/beta)*(alpha/beta)^alpha*A;
          -(alpha*exp(-alpha)*(alpha/beta)^alpha*A)/beta];
    
end
end