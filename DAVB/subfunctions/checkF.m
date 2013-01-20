function F = checkF(posterior,out)
options = out.options;
priors = options.priors;
y = out.y;
E = posterior.a_sigma./posterior.b_sigma;
V = posterior.a_sigma./posterior.b_sigma^2;
E0 = priors.a_sigma./priors.b_sigma;
V0 = priors.a_sigma./priors.b_sigma^2;
ElogS = psi(posterior.a_sigma) - log(posterior.b_sigma);
P = posterior.muPhi;
dP = P-priors.muPhi;
np = length(find(diag(priors.SigmaPhi)~=0));
F = - 0.5*np*log(2*pi) ...
    - 0.5*VBA_logDet(priors.SigmaPhi) ...
    - 0.5*dP'*VB_inv(priors.SigmaPhi)*dP ...
    - VB_KL(E,V,E0,V0,'Gamma') ...
    + 0.5*VBA_logDet(posterior.SigmaPhi) ...
    + 0.5*np*log(2*pi);
for t=1:out.dim.n_t
    ny = length(find(diag(priors.iQy{t})~=0));
    dy = y(:,t)-out.suffStat.gx(:,t);
    F = F ...
        - 0.5*ny*log(2*pi) ...
        + 0.5*VBA_logDet(priors.iQy{t}) ...
        + 0.5*ny*ElogS ...
        - 0.5*E*(dy'*priors.iQy{t}*dy); 
end
% 
% SSE = E*(dy'*priors.iQy{t}*dy)
% dF = + 0.5*ny*ElogS-VB_KL(E,V,E0,V0,'Gamma')
% S = 0.5*VBA_logDet(posterior.SigmaPhi) + 0.5*np*log(2*pi)
% ldQ = VBA_logDet(options.priors.iQy{t}) - VBA_logDet(priors.SigmaPhi)
