function [ehat,v_e,etahat,v_eta] = VBA_getNoise(posterior,out)
% returns the Laplace approximation to the innovations posterior density
% function [ehat,v_e,etahat,v_eta] = VBA_getNoise(posterior)
% This may be useful if one is interested in recovering, e.g., the state
% noise that enters and perturbs the system. In particular, posterior
% covariances are attached to state noise, which means that statistical
% inference can be performed in the usual way...
% IN:
%   - posterior/out: output structures of VBA_NLStateSpaceModel.m
% OUT:
%   - ehat/etahat: the 1st-order moment of the posterior density on the
%   measurement (resp. state) noise
%   - v_e/v_eta: the second-order moment of the posterior density on the
%   measurement (resp. state) noise

dim = out.dim;
options = out.options;
u = out.u;
y = out.y;

ehat = zeros(dim.p,dim.n_t);
v_e = cell(1,dim.n_t);
etahat = zeros(dim.n,dim.n_t);
v_eta = cell(1,dim.n_t);

% initial condition
if dim.n > 0
    [fx,dfdx,dfdp] = VBA_evalFun('f',posterior.muX0,posterior.muTheta,u(:,1),options,dim,1);
    etahat(:,1) = posterior.muX(:,1) - fx;
    if isinf(posterior.a_alpha) && isequal(posterior.b_alpha,0)
        v_eta{1} = zeros(dim.n,dim.n);
    else
        v_eta{1} = posterior.SigmaX.current{1} ...
            + dfdx'*posterior.SigmaX0*dfdx;
        if dim.n_theta > 0
            v_eta{1} = v_eta{1} + dfdp'*posterior.SigmaTheta*dfdp;
        end
    end
end
[gx,dgdx,dgdp] = VBA_evalFun('g',posterior.muX(:,1),posterior.muPhi,u(:,1),options,dim,1);
ehat(:,1) = y(:,1) - gx;
v_e{1} = zeros(dim.p,dim.p);
if dim.n > 0
    v_e{1} = dgdx'*posterior.SigmaX.current{1}*dgdx;
end
if dim.n_phi > 0
    v_e{1} = v_e{1} + dgdp'*posterior.SigmaPhi*dgdp;
end

% loop over time samples
for t = 2:dim.n_t
    if dim.n > 0
        [fx,dfdx,dfdp] = VBA_evalFun('f',posterior.muX(:,t-1),posterior.muTheta,u(:,t),options,dim,t);
        etahat(:,t) = posterior.muX(:,t) - fx;
        if isinf(posterior.a_alpha) && isequal(posterior.b_alpha,0)
            v_eta{t} = zeros(dim.n,dim.n);
        else
            P = [-dfdx',eye(dim.n)];
            jointCov = ...
                [ posterior.SigmaX.current{t}  posterior.SigmaX.inter{t-1}'
                posterior.SigmaX.inter{t-1}  posterior.SigmaX.current{t-1} ];
            v_eta{t} = P*jointCov*P';
            if dim.n_theta > 0
                v_eta{t} = v_eta{t} + dfdp'*posterior.SigmaTheta*dfdp;
            end
        end
    end
    [gx,dgdx,dgdp] = VBA_evalFun('g',posterior.muX(:,t),posterior.muPhi,u(:,t),options,dim,t);
    ehat(:,t) = y(:,t) - gx;
    v_e{t} = zeros(dim.p,dim.p);
    if dim.n > 0
        v_e{t} = dgdx'*posterior.SigmaX.current{t}*dgdx;
    end
    if dim.n_phi > 0
        v_e{t} = v_e{t} + dgdp'*posterior.SigmaPhi*dgdp;
    end
end



