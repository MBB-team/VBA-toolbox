function [gx] = g_conv_approx(~,P,u,in)

gx = P(end); % constant as starting point
nt = round(size(u,1)./in.dim.nu);

% dgdx = [];

% persistent K;
% if size(K,1) ~= nt
%     K = zeros(nt,nt);
% end

% dgdp = zeros(size(P,1),nt);



for i=1:in.dim.nu 
    Pi = P((i-1)*3+1:3*i);
    kernel = kernel_sinexp(Pi(1),Pi(2),Pi(3),0:in.deltat:in.dim.n_t);
    ui = u((i-1)*nt+1:i*nt);
    ki = conv(ui,kernel);
    gx = gx + ki(1:numel(ui));
%     for j=1:nt
%         idx = j:min(j+in.dim.n_t-1,nt) ;
%         K(idx,j) = kernel(1:min(in.dim.n_t,nt-j+1)) ;
%         dgdp((i-1)*3+1:3*i,idx) = dgdp((i-1)*3+1:3*i,idx) + [ui(idx),ui(idx),ui(idx)]'.*dkdp(:,1:min(in.dim.n_t,nt-j+1)) ; 
%     end
%     gx = gx + K*ui;
end
% dgdp(end,:) = 1;

end

function [k,dkdp,landmarks,dadp]=kernel_sinexp(A,phi,beta,tList) % 


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