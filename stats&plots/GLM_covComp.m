function [theta,e,hy,hp] = GLM_covComp(y,X,Qy,Qp,my,mp)
% wraps GLM with multiple covariance components estimation into VBA
% function [theta,e] = GLM_covComp(y,X,Qy,Qp,my,mp)
% IN:
%   - y: px1 data vector
%   - X: pxn_phi GLM design matrix
%   - Qy: nqyx1 cell array of data covariance components
%   - Qp: nqpx1 cell array of GLM coef covariance components
%   - my: prior mean on data precision hyperparameters in log-space
%   - mp: prior mean on GLM coef precision hyperparameters in log-space
% OUT:
%   - theta: n_phix1 vector of estimated GLM coef
%   - e: px1 vector of estimated residuals

nqy = length(Qy);
nqp = length(Qp);

[ny,np] = size(X);

dim.p = size(y,1);
dim.n_phi = (np+1)*nqp + (ny+1)*nqy;
dim.n_theta = 0;
dim.n = 0;
dim.n_t = 1;

priors.muPhi = zeros(dim.n_phi,1);
priors.SigmaPhi = zeros(dim.n_phi,dim.n_phi);

ind = cell(4,1);
last = 0;
for i=1:nqy
    ind{1}(i) = last+1;
    ind2i = ind{1}(i)+1:ind{1}(i)+ny;
    ind{2} = [ind{2};ind2i];
    last = ind2i(end);
    priors.muPhi(ind{1}(i)) = my(i);
    priors.SigmaPhi(ind{1}(i),ind{1}(i)) = 4;
    priors.SigmaPhi(ind2i,ind2i) = Qy{i};
end
for j=1:nqp
    ind{3}(j) = last+1;
    ind4j = ind{3}(j)+1:ind{3}(j)+np;
    ind{4} = [ind{4};ind4j];
    last = ind4j(end);
    priors.muPhi(ind{3}(j)) = mp(j);
    priors.SigmaPhi(ind{3}(j),ind{3}(j)) = 4;
    priors.SigmaPhi(ind4j,ind4j) = Qp{j};
end
priors.a_sigma = 1e8;
priors.b_sigma = 1e0;

inG.X = X;
inG.ind = ind;
inG.dim = struct('ny',ny,'np',np,'nqy',nqy,'nqp',nqp);

options = struct(...
    'priors',priors,...
    'inG',inG,...
    'updateHP',0,...
    'MaxIter',128,...
    'GnMaxIter',64);
g_fname = @g_wrapGLM;

[posterior,out] = VBA_NLStateSpaceModel(y,[],[],g_fname,dim,options);


[theta,e,hy,hp] = collapseGLM(posterior.muPhi,inG);





function [gx] = g_wrapGLM(x,P,u,in)
[theta,e] = collapseGLM(P,in);
gx = in.X*theta + e;

function [theta,e,hy,hp] = collapseGLM(P,in)
P1 = P(in.ind{1});
P2 = reshape(P(in.ind{2}),in.dim.ny,in.dim.nqy);
P3 = P(in.ind{3});
P4 = reshape(P(in.ind{4}),in.dim.np,in.dim.nqp);
hp = exp(P3);
theta = P4*hp;
hy = exp(P1);
e = P2*hy;

