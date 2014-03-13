function [ddydphi,d2gdx2,logL,dy,dy2,vy]=VBA_get_dL(gx,dG_dPhi,y,type,Qy,sigmaHat)

if nargin<5
end

% prediction error
dy = y - gx;

switch type
    
    case 0 %--- normal
        vy=(1./sigmaHat).*diag(VBA_inv(Qy,[]));
        ddydphi = sigmaHat.*(dG_dPhi*Qy*dy);
        d2gdx2 = sigmaHat.*(dG_dPhi*Qy*dG_dPhi');
        dy2=dy'*Qy*dy ;
        logL= -0.5*sigmaHat.*dy2;
        
    case 1 %--- binomial
        gx = checkGX_binomial(gx);
        vy = gx.*(1-gx) ;
        ddydphi = dG_dPhi*(dy./vy);
        temp = y./(gx).^2 - (y-1)./(1-gx).^2;
        d2gdx2 = dG_dPhi*diag(temp)*dG_dPhi';    
        logL = y'*log(gx) + (1-y)'*log(1-gx);
        dy2 = sum(temp);
        
   case 2   %--- multinomial       
        gx = checkGX_binomial(gx);
        vy = gx.*(1-gx) ;
        ddydphi = dG_dPhi*(y./gx);
        d2gdx2 = dG_dPhi*diag(y./gx.^2)*dG_dPhi';
        dy2 = sum(y./(gx).^2);
        logL = log(gx)'*y;
   
end



end