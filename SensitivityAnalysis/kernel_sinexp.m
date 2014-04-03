function [k,dkdp,landmarks,dadp]=kernel_sinexp(A,phi,beta,tList)


timeline = phi*pi*tList;
k = A*log(1+beta)*sin(timeline).*exp(-beta*timeline);


dkdp = [ log(1+beta)*sin(timeline).*exp(-beta*timeline) ;
         k.*pi.*tList.*(1./tan(timeline+eps)-beta);
         k.*(1./((beta+1)*log(beta+1))- timeline) ];

if nargout > 2
    % time to peak: alpha/beta
%     landmarks.tMax = .5/phi ;
%     landmarks.tMax=acos( sqrt(1/(1+(phi*pi/beta)^2)) )/(phi*pi);
    landmarks.tMax=atan(phi*pi/beta)/(phi*pi);
    landmarks.aMax = beta*A*exp(-beta*landmarks.tMax);
    
    dadp=[(phi*pi*exp(-(beta*atan((phi*pi)/beta))/(phi*pi)))/sqrt((phi^2*pi^2)/beta^2+1);
          (beta*abs(beta)*exp(-(beta*atan((phi*pi)/beta))/(phi*pi))*atan((phi*pi)/beta)*A)/(phi*sqrt(phi^2*pi^2+beta^2));
          -(abs(beta)*exp(-(beta*atan((phi*pi)/beta))/(phi*pi))*(beta*atan((phi*pi)/beta)-phi*pi)*A)/(beta*sqrt(phi^2*pi^2+beta^2))];
    
end
end