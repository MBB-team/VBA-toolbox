function [suffStat] = VBA_getVarEnergy(suffStat,posterior)


% Get precision parameters
sigmaHat = posterior.a_sigma./posterior.b_sigma;
alphaHat = posterior.a_alpha./posterior.b_alpha;

% observation parameters
indIn = options.params2update.phi;
Q = options.priors.SigmaPhi(indIn,indIn);
iQ = VB_inv(Q,[]);
dP = posterior.muPhi-options.priors.muPhi;
Iphi = -0.5.*dP(indIn)'*iQ*dP(indIn) -0.5*sigmaHat.*suffStat.dy2;
if ~options.ignoreMF
    Iphi = Iphi -0.5*sigmaHat*suffStat.SXd2gdx2;
end
if isweird(Iphi)
    Iphi = -Inf;
end

% evolution parameters
indIn = options.params2update.theta;
Q = options.priors.SigmaTheta(indIn,indIn);
iQ = VB_inv(Q,[]);
dP = posterior.muTheta-options.priors.muTheta;
Itheta= -0.5.*dP(indIn)'*iQ*dP(indIn) -0.5*alphaHat.*suffStat.dx2;
if ~options.ignoreMF
    Itheta = Itheta -0.5*alphaHat*suffStat.SXd2fdx2;
end
if isweird(Itheta)
    Itheta = -Inf;
end

% hidden states
IX = -0.5.*sigmaHat.*suffStat.dy2 -0.5*alphaHat.*suffStat.dx2;

