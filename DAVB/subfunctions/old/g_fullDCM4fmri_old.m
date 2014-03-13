function [gx,dG_dX,dG_dPhi] = g_fullDCM4fmri(Xt,Phi,ut,inG)
% [broken]

% function [gx,dG_dX,dG_dPhi] = g_fullDCM4fmri(Xt,Phi,ut,inG)
% This function evaluates the observation function of the neuronal level in
% DCM for fMRI. Note that this "generalized" observation function includes
% the HRF convolution model (balloon model).


persistent xh dxhdx dxhdp ph ti
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

% Check whether the system is at initial state
if isempty(xh)
    ti = 0;
    for i=1:inG.n
        xh{i} = [0;1;1;1];
        dxhdx{i} = zeros(1,4);
        dxhdp{i} = zeros(length(inG.indHRF{i}),4);
    end
end

if isempty(Phi)    %- Function is called from f_fullDCM4fmri
    
    % Get current observation parameters mode
    if ~isempty(ph)
        Phi = ph;
    else % only the first iteration
        Phi = inG.Phi;
    end

    %- hemodynamic response convolution operation
    n = size(Xt,1);
    for i=1:n
        % hemodynamic states evolution
        P = Phi(inG.indHRF{i});
        [fx,dfdxh,dfdp] = f_HRF(xh{i},P,Xt(i),inG);
        xh{i} = fx;
        % update gradients
        dxhdx{i} = dxhdx{i}*dfdxh + [1 0 0 0];
        dxhdp{i} = dfdp + dxhdp{i}*dfdxh;
    end

else    %- Function is called outside of f_fullDCM4fmri
    
    % Update persistent observation parameters
    ph = Phi;
    ti = ti+1;
    
    %- Get predicted observations and gradients
    n = size(Xt,1);
    gx = zeros(n,1);
    dG_dX = zeros(n,n);
    dG_dPhi = zeros(size(Phi,1),n);
    for i=1:n
        % static observation function
        P = Phi(inG.indHRF{i});
        [gx(i),dgdxh,dgdp] = g_HRF(xh{i},P,ut,inG);
        % get gradients
        dG_dX(i,:) = dxhdx{i}*dgdxh;
        dG_dPhi(inG.indHRF{i},i) = dgdp + dxhdp{i}*dgdxh;
%         % empty gradients w.r.t. states
%         dxhdx{i} = zeros(1,4);
    end

end

if ti==inG.n_t
    xh = [];
end



