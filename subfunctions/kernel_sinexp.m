function [k,dkdp,landmarks,dadp]=kernel_sinexp(A,phi,beta,tList)


base = beta*sin(phi*pi*tList).*exp(-beta*tList);
k=A.*base;


dkdp = [ base;
         beta*A*pi*tList.*exp(-beta*tList).*cos(phi*pi*tList);
         -(beta*tList-1).*exp(-beta*tList).*sin(phi*pi*tList)*A ];

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