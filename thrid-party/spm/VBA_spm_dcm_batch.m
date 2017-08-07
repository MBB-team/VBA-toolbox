function [err,F,DCMnames] = VBA_spm_dcm_batch(DCM0dir,newVOIdir,newSPMdir,Unames,Ynames,est)
% replicates one single subject DCM analysis over a set of subjects
% function [err,F] = spm_dcm_batch(DCM0dir,newVOIdir,newSPMdir,Unames,Ynames)
% IN:
%   - DCM0dir: directory where the original DCM files to be cloned are
%   stored
%   - newVOIdir: nsX1 cell-array of directories where the VOI files are
%   stored
%   - newSPMdir: nsX1 cell-array of directories where the SPM files are
%   stored
%   - Unames: nsX1 cell array of ndcmX1 cell arrays of input names. For
%   each of the ns subjects, the function will look for each one of them in
%   the SPM file. NB: the order should match the original DCM structure. If
%   empty or unspecified, these will be set to the names of the inputs in
%   the original DCM structure
%   - Ynames: nsX1 cell array of ndcmX1 cell arrays of ROI names. For
%   each of the ns subjects, the function will look for each one of them in
%   the VOI files. NB: the order should match the original DCM structure.
%   If empty or unspecified, these will be set to the names of the ROIs in
%   the original DCM structure.
%   - est: flag for DCM estimation (1: estimate, 0: do not estimate)
% OUT:
%   - err: nsXndcm array of error flags. It is 0 if DCM appropriately
%   cloned, 1 problem with ROIs, 2 if problem with inputs, 3 if problem
%   with the original DCM and 4 if problem when saving the DCM structure, 5
%   if DCM could not be inverted and 6 if it could not be saved.
%   - F: nsXndcmX2 array of free energies.
%   - DCMnames: ndcmX1 cell array of DCM names

ns = length(newVOIdir);
if ~exist('Unames','var') || isempty(Unames)
    for i=1:ns
        Unames{i} = [];
    end
end
if ~exist('Ynames','var') || isempty(Ynames)
    for i=1:ns
        Ynames{i} = [];
    end
end
if ~exist('est','var') || isempty(est)
    est = 0;
end

% 1- get original DCM file names
[DCM0fnames] = getFilesInDir(DCM0dir,'DCM_');

% 2- loop over new subjects and clone the original DCMs
% NB: this assumes that all DCMs have identical inputs/VOIs
ndcm = length(DCM0fnames);

err = zeros(ns,ndcm);
F = zeros(ns,ndcm,2);
DCMnames = cell(ndcm,1);
for i=1:ns
    disp(['Looking in ',newVOIdir{i},' ...'])
    % 2.1 get VOI file names
    [VOIfnames] = getFilesInDir(newVOIdir{i},'VOI_');
    [SPMfname] = getFilesInDir(newSPMdir{i},'SPM');
    % 2.2 clone each original DCM
    for j=1:ndcm
        lo = load(DCM0fnames{j});
        [pathstr,name,ext] = fileparts(DCM0fnames{j});
        DCMfname = [newVOIdir{i},filesep,name,'.mat'];
        if i==1
            DCMnames{j} = name;
        end
        skip  = 0;
        if exist(DCMfname,'file')
            lo = load(DCMfname);
            if isfield(lo.DCM,'VBA')
                skip = 1;
                disp(['    - skipping ',name,' ...'])
                try
                    F(i,j,1) = lo.DCM.VBA.out.F;
                    F(i,j,2) = lo.DCM.VBA.out.options.init.out.F;
                catch
                    F(i,j,1) = NaN;
                    F(i,j,2) = NaN;
                end
            end
        end
        if ~skip
            disp(['    - cloning ',name,' ...'])
            [DCM,err(i,j)] = VBA_spm_dcm_clone(lo.DCM,VOIfnames,SPMfname{1},DCMfname,Unames{i},Ynames{i});
            if ~err(i,j) && est
                try
                    disp(['    - inverting ',name,' ...'])
                    DCM.VBA = 'idle';
                    if VBA_spm_check_version('matlab','7') >= 0
                        save(DCMfname,'-V6','DCM');
                    else
                        save(DCMfname,'DCM');
                    end
                    [DCM] = VBA_spm_sdcm_estimate(DCM,1,0,0);
                    F(i,j,1) = DCM.VBA.out.F;
                    F(i,j,2) = DCM.VBA.out.options.init.out.F;
                catch
                    disp('       Error: could not invert DCM')
                    err(i,j) = 5;
                end
                try
                    disp(['    - saving ',name,' ...'])
                    if VBA_spm_check_version('matlab','7') >= 0
                        save(DCMfname,'-V6','DCM');
                    else
                        save(DCMfname,'DCM');
                    end
                catch
                    disp('       Error: could not save DCM')
                    err(i,j) = 6;
                end
            end
            
        end
    end
end


function [flist] = getFilesInDir(dirName,prefix)
% list files in dirName, whose name begins with prefix
fn = dir(dirName);
flist = cell(0);
for i=1:length(fn)
    if length(fn(i).name) >= length(prefix)
        tmp = fn(i).name(1:length(prefix));
        if isequal(tmp,prefix)
            flist{end+1} = [dirName,filesep,fn(i).name];
        end
    end
end


