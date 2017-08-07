function [sF,sDCM] = VBA_spm_split_dcm(DCM,nsplits)
% split DCM inversion using time series partition
% function [sF,sDCM] = spm_split_dcm(DCM,nsplits)
% IN:
%   - DCM: the original DCM structure
%   - nsplits: the number of splits in the partition
% OUT:
%   - sF: log-evidences for each split in the partition
%   - sDCM: cell array of DCM structures for each split

try; nsplits; catch; nsplits=4; end

% 1- split DCM structure
[sDCM] = split_DCM(DCM,nsplits);

% 2- invert split-models
sF = zeros(nsplits,1);
for i=1:nsplits
    try
        % fprintf(1,'confound removal...')
        nq = size(sDCM{i}.Y.X0,2);
        sDCM{i}.Y.X0 = sDCM{i}.Y.X0(ceil(nq./nsplits),:);
        sDCM{i}.Y.X0 = ones(sDCM{i}.v,1);
        % sDCM{i}.delays =0.*sDCM{i}.delays;
        [sDCM{i}] = VBA_spm_dcm_estimate(sDCM{i});
        sF(i) = sDCM{i}.F;
        fprintf(1,' OK.\n')
    catch
        sF(i) = NaN;
    end
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

