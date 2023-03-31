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

X = [posterior.muX0 posterior.muX];

%% loop over time samples

for t = 1 : dim.n_t

    % evolution noie
    % =====================================================================
    if dim.n > 0
        
        % (1st moment)
        % -----------------------------------------------------------------
        [fx,dfdx,dfdp] = VBA_evalFun('f',X(:,t),posterior.muTheta,u(:,t),options,dim,t);
        etahat(:,t) = X(:,t+1) - fx;
        
        % (2nd moment)
        % -----------------------------------------------------------------
        % deterministic system
        if isinf(posterior.a_alpha) && isequal(posterior.b_alpha,0)     
            v_eta{t} = zeros(dim.n,dim.n);

        % stochastic system
        else
            if t == 1
                v_eta{t} = posterior.SigmaX.current{t} ...
                    + dfdx'*posterior.SigmaX0*dfdx;
            else
                P = [-dfdx',eye(dim.n)];
                jointCov = ...
                    [ posterior.SigmaX.current{t}  posterior.SigmaX.inter{t-1}'
                    posterior.SigmaX.inter{t-1}  posterior.SigmaX.current{t-1} ];
                v_eta{t} = P*jointCov*P';
            end
            if dim.n_theta > 0
                v_eta{t} = v_eta{t} + dfdp'*posterior.SigmaTheta*dfdp;
            end
        end
    end

    % observation
    % =====================================================================

    % (1st moment)
    % ---------------------------------------------------------------------
    [gx,dgdx,dgdp] = VBA_evalFun('g',X(:,t+1),posterior.muPhi,u(:,t),options,dim,t);
    ehat(:,t) = y(:,t) - gx;

    % (2nd moment)
    % ---------------------------------------------------------------------
    v_e{t} = zeros(dim.p,dim.p);
    if dim.n > 0
        v_e{t} = dgdx'*posterior.SigmaX.current{t}*dgdx;
    end
    if dim.n_phi > 0
        v_e{t} = v_e{t} + dgdp'*posterior.SigmaPhi*dgdp;
    end
end



