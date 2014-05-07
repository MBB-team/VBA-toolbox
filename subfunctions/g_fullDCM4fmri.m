function [gx,dG_dX,dG_dPhi] = g_fullDCM4fmri(Xt,Phi,ut,inG)
% full DCM observation function (embedding dynamic Ballon model)

% function [gx,dG_dX,dG_dPhi] = g_fullDCM4fmri(Xt,Phi,ut,inG)
% This function evaluates the observation function of the neuronal level in
% DCM for fMRI. Note that this "generalized" observation function includes
% the HRF convolution model (balloon model).


persistent xh dxhdx dxhdp ph ti ko
% These variables are described bellow:
%   - xh{i}: current hemodynamic states for region i.
%   - dxhdx{i}: derivatives of the hemodynamic states of region i w.r.t.
%   neuronal state of region i.
%   - dxhdp{i}: derivatives of the hemodynamic states of region i w.r.t.
%   hemodynamic parameters of region i.
%   - ph: the last remembered hemodynamic parameters. This is only used
%   when this function is called from f_fullDCM4fmri.m. It is updated each
%   time this function is called outside f_fullDCM4fmri.m.
%   - ti: the index within the time series. This is used to reinitialize
%   the persistent variables.
% These are reinitialized at the beginning of the time series.

try % just to get lastly updated phi
    
    inG.getPhi;
    ph = Phi;
    
catch

    n = size(Xt,1);

    % Check whether the system is at initial state
    if isempty(xh)
        ti = 0;
        xh = kron(ones(n,1),[0;0;0;0]);
        dxhdx = zeros(n,n*4);
        dxhdp = zeros(size(inG.Phi,1),n*4);
        ko = kron(eye(n),[1 0 0 0]);
    end


    if isempty(Phi)    %- Function is called from f_fullDCM4fmri

        % Get current observation parameters mode
        if ~isempty(ph)
            Phi = ph;
        else % only the first iteration
            Phi = inG.Phi;
        end

        %- hemodynamic response convolution operation
        [fx,dfdxh,dfdp] = f_HRF2(xh,Phi,Xt,inG);

        xh = fx;
        % update gradients
        dxhdx = dxhdx*dfdxh + ko;
        dxhdp = dfdp + dxhdp*dfdxh;

    else    %- Function is called outside of f_fullDCM4fmri

        % Update persistent observation parameters
        ph = Phi;
        ti = ti+1;

        %- Get predicted observations and gradients
        [gx,dgdxh,dgdp] = g_HRF3(xh,Phi,ut,inG);
        %
        %     dgdxh = numericDiff('g_HRF',1,xh,Phi,ut,inG);
        %     dgdp = numericDiff('g_HRF',2,xh,Phi,ut,inG);

        dG_dX = dxhdx*dgdxh;
        dG_dPhi = dgdp + dxhdp*dgdxh;
        %         % empty gradients w.r.t. states
        %         dxhdx = zeros(n,n*4);

    end

    if ti==inG.n_t
        xh = [];
    end
end


