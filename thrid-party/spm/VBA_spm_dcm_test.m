function DCM = VBA_spm_dcm_test(DCM)
% tests whether DCM parameters are zero using Savage-Dickey ratios
% function DCM = spm_dcm_test(DCM)
% IN:
%   - DCM: the DCM structure
% OUT:
%   - DCM: the DCM structure, augmented with field .sdr, which contains the
%   Savage-Dickey ratios in the same format as DCM.Ep

n = size(DCM.Ep.A,1);
nu = size(DCM.Ep.C,2);

E = VBA_spm_vec(DCM.Ep);
V = DCM.Cp;
E0 = VBA_spm_vec(DCM.M.pE);
V0 = DCM.M.pC;

sdr.A = zeros(n,n);
for i=1:n
    for j=1:n
        % 0- identify param index in DCM.Cp
        tmp = DCM.Vp;
        tmp.A(i,j) = 0;
        ind = find(full(VBA_spm_vec(tmp))~=full(VBA_spm_vec(DCM.Vp)));
        % 1- define reduced model
        E0r = E0;
        E0r(ind) = 0;
        V0r = V0;
        V0r(ind,:) = 0;
        V0r(:,ind) = 0;
        % 2- call Savage-Dickey ratios
        if isempty(ind) % this param was fixed
            if isequal(E0,E0r)
                dF = Inf;
            else
                dF = -Inf;
            end
        else
            dF = VBA_spm_log_evidence(E,V,E0,V0,E0r,V0r);
        end
        % 3- store Bayes factor
        sdr.A(i,j) = 1./(1+exp(dF));
    end
end

sdr.B = zeros(n,n,nu);
for i=1:n
    for j=1:n
        for k=1:nu
            tmp = DCM.Vp;
            tmp.B(i,j,k) = 0;
            ind = find(full(VBA_spm_vec(tmp))~=full(VBA_spm_vec(DCM.Vp)));
            E0r = E0;
            E0r(ind) = 0;
            V0r = V0;
            V0r(ind,:) = 0;
            V0r(:,ind) = 0;
            if isempty(ind) % this param was fixed
                if isequal(E0,E0r)
                    dF = Inf;
                else
                    dF = -Inf;
                end
            else
                dF = VBA_spm_log_evidence(E,V,E0,V0,E0r,V0r);
            end
            sdr.B(i,j,k) = 1./(1+exp(dF));
        end
    end
end


sdr.C = zeros(n,nu);
for i=1:n
    for j=1:nu
        tmp = DCM.Vp;
        tmp.C(i,j) = 0;
        ind = find(full(VBA_spm_vec(tmp))~=full(VBA_spm_vec(DCM.Vp)));
        E0r = E0;
        E0r(ind) = 0;
        V0r = V0;
        V0r(ind,:) = 0;
        V0r(:,ind) = 0;
        if isempty(ind) % this param was fixed
            if isequal(E0,E0r)
                dF = Inf;
            else
                dF = -Inf;
            end
        else
            dF = VBA_spm_log_evidence(E,V,E0,V0,E0r,V0r);
        end
        sdr.C(i,j) = 1./(1+exp(dF));
    end
end

nD = size(DCM.Ep.D,3);
if nD > 0
    for i=1:n
        for j=1:n
            for k=1:n
                tmp = DCM.Vp;
                tmp.D(i,j,k) = 0;
                ind = find(full(VBA_spm_vec(tmp))~=full(VBA_spm_vec(DCM.Vp)));
                E0r = E0;
                E0r(ind) = 0;
                V0r = V0;
                V0r(ind,:) = 0;
                V0r(:,ind) = 0;
                if isempty(ind) % this param was fixed
                    if isequal(E0,E0r)
                        dF = Inf;
                    else
                        dF = -Inf;
                    end
                else
                    dF = VBA_spm_log_evidence(E,V,E0,V0,E0r,V0r);
                end
                sdr.D(i,j,k) = 1./(1+exp(dF));
            end
        end
    end
end


DCM.sdr = sdr;

