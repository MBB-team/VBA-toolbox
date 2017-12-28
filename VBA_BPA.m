function [p_BPA] = VBA_BPA(priors0,posteriors0)
% performs Bayesian Parameters Averaging (BPA)
% function [p,o] = VBA_BPA(posterior,F)
% IN:
%   - priors0: a Kx1 cell-array of VBA prior structures (where K is the number of subjects)
%   - posteriors0: a Kx1 cell-array of VBA posterior structures (where K is the number of subjects)
% OUT:
%   - p_BPA: the resulting posterior structure, with the first two moments of
%   the group-level probability density functions

p0 = priors0;
p1 = posteriors0;

% observation parameters
[p_BPA.muPhi,p_BPA.SigmaPhi] = get2moments(p0,p1,'Phi');

% evolution parameters
[p_BPA.muTheta,p_BPA.SigmaTheta] = get2moments(p0,p1,'Theta');

% initial conditions
[p_BPA.muX0,p_BPA.SigmaX0] = get2moments(p0,p1,'X0');

% hidden states
p_BPA.muX=[];
p_BPA.SigmaX=[];

% data precision
p_BPA.b_sigma =[];
p_BPA.a_sigma =[];

% hidden state precision
p_BPA.b_alpha = [];
p_BPA.a_alpha =[];

end

function [muG,sigmaG] = get2moments(p0,p1,paramType)

    K = length(p0); % # subjects


    % define group priors 
    mu0 = p0{1}.(['mu' paramType]);
    sigma0 = p0{1}.(['Sigma' paramType]);
    muG = mu0;
    sigmaG = sigma0;
    muSub = zeros(numel(mu0),K);
    sigmaSub =  cell(1,K);
    a0 = ones(numel(mu0),1);
    b0 = ones(numel(mu0),1);
    aG = a0;
    bG = b0;
    ind_ffx = find(infLimit(a0,b0)==1);
    ind_in = find(diag(sigma0)~=0);
    
    % loop across subjects
    for k=1:K

        % subject-level posterior
        mu = p1{k}.(['mu' paramType]);
        sigma = p1{k}.(['Sigma' paramType]);
        
%         % update
%         tempSigma = inv( inv(sigmaG) + inv(sigma) );
%         muG = tempSigma*inv(sigmaG)*muG + tempSigma*inv(sigma)*mu;
%         sigmaG = tempSigma ;
        
        % store
        muSub(:,k) = mu;
        sigmaSub{k} =  sigma;
        
    end
    
    % VB-updating
    [muG,sigmaG,aG,bG] = MFX_VBupdate(muG,VBA_inv(sigmaG),...
                                 muSub,sigmaSub,...
                                 aG,bG,...
                                 a0,b0,...
                                 ind_ffx,ind_in);
    


end

function [m,V,a,b] = MFX_VBupdate(m0,iV0,ms,Vs,a,b,a0,b0,indffx,indIn)
ns = size(ms,2);
n = size(m0,1);
sm = 0;
sv = 0;
wsm = 0;
sP = 0;
indrfx = setdiff(1:n,indffx);
indrfx = intersect(indrfx,indIn);
indffx = intersect(indffx,indIn);
iQ = diag(a(indrfx)./b(indrfx));
for i=1:ns
    % RFX
    sm = sm + ms(indrfx,i);
    e = ms(indrfx,i)-m0(indrfx);
    sv = sv + e.^2 + diag(Vs{i}(indrfx,indrfx));
    % FFX
    tmp = VBA_inv(Vs{i});
    wsm = wsm + tmp*ms(:,i);
    sP = sP + tmp;
end
% RFX
V = zeros(n,n);
m = m0;
V(indrfx,indrfx) = VBA_inv(iV0(indrfx,indrfx)+ns*iQ);
m(indrfx) = V(indrfx,indrfx)*(iV0(indrfx,indrfx)*m0(indrfx)+iQ*sm);
a(indrfx) = a0(indrfx) + 0.5*ns;
b(indrfx) = b0(indrfx) + 0.5*(sv(indrfx)+ns*diag(V(indrfx,indrfx)));
% FFX
if ~isempty(indffx)
    tmp = VBA_inv(sP);
    V(indffx,indffx) = tmp(indffx,indffx);
    m(indffx) = V(indffx,indffx)*wsm(indffx);
end
end

function il = infLimit(a,b)
il = isinf(a).*isequal(b,0);
end






