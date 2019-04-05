function [ddydphi, d2gdx2, logL, dy, dy2, vy]= VBA_get_dL (gx, dG_dPhi, y, type, Qy, sigmaHat)
% // VBA toolbox //////////////////////////////////////////////////////////
%
% [ddydphi, d2gdx2, logL, dy, dy2, vy]= VBA_get_dL (gx, dG_dPhi, y, type, Qy, sigmaHat)
% Compute usefull intermediate values describing the misfit between a model
% prediction and an observation.
%
% IN:
%   - gx: model prediction about the observation y (1st order moment)
%   - dG_dPhi: derivative of gx wrt observation parameters
%   - y: actual observation
%   - type: flag defining the distribution of the observation (0: gaussian,
%     1: Bernouilli, 2: categorical)
%   - Qy: scaling matrix of the 2nd order moment of the prediction (only 
%     gaussian observations)
%   - sigmaHat: scaling factor of the 2nd order moment of the prediction 
%     (only gaussian observations)
%
% OUT:
%   - ddydphi: gradient of the prediction error wrt observation parameters
%   - d2gdx2: hessian of the prediction
%   - logL: log-likelihood of the observation given the prediction
%   - dy: prediction error
%   - dy2: normalized squared deviation 
%   - vy: prediction variance
%     
% /////////////////////////////////////////////////////////////////////////

if nargin<5
end

% prediction error
dy = y - gx;

switch type
    
    case 0 %--- normal
        vy=(1./sigmaHat).*diag(VBA_inv(Qy));
        ddydphi = sigmaHat.*(dG_dPhi*Qy*dy);
        d2gdx2 = sigmaHat.*(dG_dPhi*Qy*dG_dPhi');
        dy2=dy'*Qy*dy ;
        logL = - 0.5*sigmaHat.*dy2 ;
        logL = logL + 0.5*VBA_logDet(Qy*sigmaHat) - 0.5*numel(dy2)*log(2*pi) ;
        
    case 1 %--- binomial
        gx = VBA_finiteBinomial (gx);
        vy = gx.*(1-gx) ;
        ddydphi = dG_dPhi*(dy./vy);
        temp = y./(gx).^2 - (y-1)./(1-gx).^2;
        d2gdx2 = dG_dPhi*diag(temp)*dG_dPhi';    
        logL = y'*log(gx) + (1-y)'*log(1-gx);
        dy2 = sum(temp);
        
   case 2   %--- multinomial       
        gx = VBA_finiteBinomial (gx);
        vy = gx.*(1-gx) ;
        ddydphi = dG_dPhi*(y./gx);
        d2gdx2 = dG_dPhi*diag(y./gx.^2)*dG_dPhi';
        dy2 = sum(y./(gx).^2);
        logL = log(gx)'*y;
   
end



end