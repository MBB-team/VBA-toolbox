function [gx,dgdx,dgdp] = g_HRF_distributed(Xt,P,ut,in)
% observation function for distributed HRF Balloon model (DCM for fMRI)
% function [gx,dgdx,dgdp] = g_HRF(Xt,P,ut,in)
% This function evaluates the hemodynamic static observation equation
% function.

if isfield(in,'fullDCM') && in.fullDCM
    ind_profile = in.ind3;
    ind_hrf = 1:size(P,1);
else
    ind_profile = in.ind_profile;
    ind_hrf = in.ind_hrf;
end

% call Ballon model
[w,dwdx,dwdp] = g_HRF3(Xt,P(ind_hrf),ut,in);

% get spatial profile parameters and mixing matrix
Phi = zeros(in.n_phi,in.n_reg);
for i=1:in.n_reg
    Phi(:,i) = P(ind_profile{i});
end
A = in.B*Phi;

% apply mixing matrix to observation function and jacobian
gx = A*w;
dgdx = dwdx*A';

% recollect gradient wrt observation parameters
dgdp = zeros(size(P,1),size(gx,1));
dgdp(ind_hrf,:) = dwdp*A';
for i=1:in.n_reg
    dgdp(ind_profile{i},:) = w(i).*in.B';
end
