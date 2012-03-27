function [b] = spm_dcm_risk_split(dcm_fnames,fam,P,nsplits)

nm = length(dcm_fnames);
sdcms = cell(nsplits,nm);
for i=1:nm
    if isstruct(dcm_fnames{i})
        DCM = dcm_fnames{i};
    else
        load(dcm_fnames{i})
    end
    [sDCM] = split_DCM(DCM,nsplits);
    sdcms(:,i) = sDCM;
end

b = zeros(nsplits,1);
for i=1:nsplits
    [DJS,b(i)] = dcm_risk(sdcms(i,:),ones(nm,1)./nm,0,fam,P);
end


function [splitDCM] = split_DCM(DCM,nsplits)

if nsplits==1
    splitDCM{1} = DCM;
    return
end

splitDCM = cell(nsplits,1);

nt = DCM.v;
subnt = floor(nt./nsplits);

ntu = size(DCM.U.u,1);
subntu = floor(ntu./nsplits);

for i=1:nsplits
    splitDCM{i} = DCM;
    in = subnt*(i-1)+1:subnt*i;
    splitDCM{i}.v = subnt;
    splitDCM{i}.Y.y = splitDCM{i}.Y.y(in,:);
    splitDCM{i}.Y.X0 = splitDCM{i}.Y.X0(in,:);
    nq = length(DCM.Y.Q);
    %     splitDCM{i}.Y.Q = eye(subnt);
    for j=1:nq
        tmp = zeros(nq*subnt,1);
        tmp(subnt*(j-1)+1:subnt*j) = 1;
        splitDCM{i}.Y.Q{j} = diag(tmp);
    end
    inu = subntu*(i-1)+1:subntu*i;
    splitDCM{i}.U.u = splitDCM{i}.U.u(inu,:);
end
