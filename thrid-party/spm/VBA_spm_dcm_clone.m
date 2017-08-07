function [DCM,err] = VBA_spm_dcm_clone(DCM0,VOIfnames,SPMfname,DCMfname,Unames,Ynames)
% clones an existing DCM structure and replaces appropriate data and inputs
% function [DCM] = VBA_spm_dcm_clone(DCM0,VOIfnames,SPMfname,DCMfname,Unames,Ynames)
% IN:
%   - DCM0: original DCM structure to be cloned
%   - VOIfnames: cell-array of ROI file names for the new DCM
%   - SPMfname: SPM file name for the new DCM
%   - DCMfname: the cloned DCM will be saved using this file name
%   - Unames: cell array of input names. The function will look for each
%   one of them in the SPM file. NB: the order should match the original
%   DCM structure. If empty or unspecified, these will be set to the names
%   of the inputs in the original DCM structure.
%   - Ynames: cell array of ROI names. The function will look for each
%   one of them in the VOI files. NB: the order should match the original
%   DCM structure. If empty or unspecified, these will be set to the names
%   of the ROIs in the original DCM structure.
% OUT:
%   - DCM: cloned DCM, with inputs and VOI relaced appropriately.
%   - err: Error flag. It is 0 if DCM appropriately cloned, 1 problem with
%   ROIs, 2 if problem with inputs, 3 if problem with the original DCM and
%   4 if problem when saving the DCM structure.

if ~exist('Unames','var') || isempty(Unames)
    Unames = DCM0.U.name;
end
if ~exist('Ynames','var') || isempty(Ynames)
    for i=1:numel(DCM0.xY)
        Ynames{i} = DCM0.xY(i).name;
    end
end

err = 0;

%==========================================================================
% Outputs
%==========================================================================
try
    % 1- find matched ROIs within the VOI files
    m = numel(VOIfnames);
    Y2find = 1:numel(Ynames);
    for i = 1:m
        p = load(VOIfnames{i},'xY');
        p.xY.VOIfile = VOIfnames{i};
        look4Y = ismember(Ynames,p.xY(1).name);
        if ~isempty(find(look4Y==1,1))
            xY(find(look4Y==1,1)) = p.xY;
            Y2find = setdiff(Y2find,find(look4Y==1,1));
        end
    end
    % 2- check ROI consistency
    if ~isempty(Y2find)
        disp('Error: could not find the following ROI:')
        for i = 1:length(Y2find)
            disp(['      - ',Ynames{Y2find(i)}]);
        end
        DCM = [];
        err = 1;
        return
    end
catch
    DCM = [];
    err = 1;
    return
end

%==========================================================================
% Inputs
%==========================================================================
try
    %-Get (nc) 'causes' or inputs U
    %----------------------------------------------------------------------
    if isequal(DCM0.U.u,zeros(DCM0.v,1)) && isequal(DCM0.U.name,{'null'})
        % endogenous input
        U.u = DCM0.U.u;
        U.name = DCM0.U.name;
    else
        % 1- find matched inputs within the SPM structure
        load(SPMfname)
        Sess   = SPM.Sess(xY(1).Sess);
        u2find = 1:size(DCM0.U.u,2); % indices of inputs in DCM0
        if ~isempty(Sess.U)
            % with stimuli
            U.dt = Sess.U(1).dt;
            u = length(Sess.U);
            U.name = {};
            U.u = [];
            for  i = 1:u
                for j = 1:length(Sess.U(i).name)
                    look4u = ismember(Unames,Sess.U(i).name{j});
                    if ~isempty(find(look4u==1,1))
                        U.u = [U.u,Sess.U(i).u(33:end,j)];
                        U.name{end+1} = Sess.U(i).name{j};
                        u2find = setdiff(u2find,find(look4u==1,1));
                    end
                end
            end
        end
        % 2- check inputs consistency
        if ~isempty(u2find)
            disp('Error: could not find the following inputs:')
            for i = 1:length(u2find)
                disp(['      - ',Unames{u2find(i)}]);
            end
            DCM = [];
            err = 2;
            return
        end
    end
catch
    DCM = [];
    err = 2;
    return
end


%==========================================================================
% Timings, model options, graph conncetions, responses
%==========================================================================
try
    %-Slice timings
    %----------------------------------------------------------------------
    delays = DCM0.delays; % assume identical delays
    
    %-Echo time (TE) of data acquisition
    %----------------------------------------------------------------------
    TE = DCM0.TE; % assume identical TE
    
    %-Model options
    %----------------------------------------------------------------------
    options = DCM0.options;
    
    %-Model structure
    %----------------------------------------------------------------------
    a = DCM0.a;
    b = DCM0.b;
    c = DCM0.c;
    d = DCM0.d;
        
    %-Response variables & confounds (NB: the data have been whitened)
    %----------------------------------------------------------------------
    n = length(xY);                      % number of regions
    v = length(xY(1).u);                 % number of time points
    Y.dt = SPM.xY.RT;
    Y.X0 = xY(1).X0;
    for i = 1:n
        Y.y(:,i) = xY(i).u;
        Y.name{i} = xY(i).name;
    end
    
    %-Precision components (one for each region): i.i.d. (because of W)
    %----------------------------------------------------------------------
    Y.Q = VBA_spm_Ce(ones(1,n)*v);
    
    %-Store all variables in DCM structure
    %----------------------------------------------------------------------
    DCM.a = a;
    DCM.b = b;
    DCM.c = c;
    DCM.d = d;
    DCM.U = U;
    DCM.Y = Y;
    DCM.xY = xY;
    DCM.v = v;
    DCM.n = n;
    DCM.TE = TE;
    DCM.delays = delays;
    DCM.options = options;
    
catch
    
    DCM = [];
    err = 3;
    return
    
end

try
    %-Save
    %----------------------------------------------------------------------
    if VBA_spm_check_version('matlab','7') >= 0
        save(DCMfname,'-V6','DCM');
    else
        save(DCMfname,'DCM');
    end
    
catch
    
    err = 4;
    
end


