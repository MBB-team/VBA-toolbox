function [out] = RFT_main(X,options,verbose)
% applies RFT to an input 1D-field X sampled on a regular lattice
% function [out] = RFT_main(X,options,verbose)
% RFT (Random Field Theory) allows one to derive the FWE-corrected p-values
% on some topological features of the field (e.g., local peaks and
% upcrossing thresholds), under the assumption of non-independance
% between neighbouring positions on the lattice.
% IN:
%   - X: Lx1 vector of sampled RF
%   - options: structure containing the set-inducing thresholds, i.e.:
%       .u: height threshold (cluster-level inference, default=2.33)
%       .k: extent threshold (set-level inference, default=6)
%       .R: sampled RF of residuals (for smoothness estimation). If empty,
%       then the smoothness is estimated directly on X.
%       .type: type of RF. Can be set to 'norm' (normal, default), 't'
%       (Student) or 'F' (Fisher).
%       .dof: degrees of freedom (only relevant for 't' or 'F' fields).
%   - verbose: a verbose flag
% OUT:
%   - out: a structure containing the following fields:
%       .peaks: summary stats of peak-level inference, i.e.:
%           .ind: npx1 vector of local maxima indices
%           .val: npx1 vector of local maxima values
%           .prft: npx1 vector of local maxima corrected p-values
%           .punc: npx1 vector of local maxima uncorrected p-values
%       where np is the number of local maxima on the field.
%       .clusters: summary stats of cluster-level inference, i.e.:
%           .ind: ncx1 cell-array of sets of indices
%           .k: ncx1 vector of clusters' spatial extent
%           .prft: ncx1 vector of clusters' corrected p-values
%       where nc is the number of upcrossing clusters on the field.
%       .set: summary stats of set-level inference, i.e.:
%           .c: number of upcrossing threshold with a size bigger than k
%           .prft: set-level p-value
%       .xc: the height threshold that corresponds to the nominal FWER
%       .Em: expected number of clusters (under H0)
%       .En: expected number of voxels per cluster (under H0)
%       .fwhm: the estimated smoothness of the field

try,options;catch,options=[];end
try,verbose;catch,verbose=0;end

OUTSTR = cell(0);
OUTSTR{1} = ['Date: ',datestr(clock)];
if verbose
    disp(' ')
    disp('-- 1D-RFT analysis --')
    disp(OUTSTR{1})  
end

if isempty(options)
    options.FWER = 0.05;
    options.u = VBA_spm_invNcdf(0.99,0,1);
    options.k = 6;
    options.R = [];
    options.type = 'norm';
    options.dof = NaN; % only relevant for 't' or 'F' fields
    options.mask = [];
    if verbose
        disp(['Using default set-inducing thresholds (X>',num2str(options.u,'%3.2f'),', k>',num2str(options.k),').'])
    end
else
    if ~isfield(options,'FWER')
        options.FWER = 0.05;
    end
    if ~isfield(options,'k')
        options.k = 6;
        if verbose
            disp(['Using default set-inducing threshold (k>',num2str(options.k),').'])
        end
    end
    if ~isfield(options,'R')
        options.R = [];
        if verbose
            disp(['Estimating smoothness on statistical map.'])
        end
    end
    if ~isfield(options,'type')
        options.type = 'norm';
    end
    if ~isfield(options,'dof')
        if ~isequal(options.type,'norm')
            disp('RFT: error: must provide valid degrees of freedom!')
            out = [];
        else
            options.dof = NaN;
        end
    end
    if ~isfield(options,'u')
        switch options.type
            case 'norm'
                options.u = VBA_spm_invNcdf(0.99,0,1);
            case 't'
                options.u = VBA_spm_invTcdf(0.99,options.dof);
            case 'F'
                options.u = VBA_spm_invFcdf(0.99,options.dof(1),options.dof(2));
        end
        if verbose
            disp(['Using default cluster-inducing threshold (X>',num2str(options.u,'%3.2f'),').'])
        end
    end
    if ~isfield(options,'mask')
        options.mask = [];
    end
end

if ~isempty(options.mask)
    if length(options.mask)~=length(X) || ~isequal(VBA_vec(unique(options.mask)),[0;1])
        disp(['RFT-1D: error: invalid mask provided!'])
        out = [];
        return
    end
end

out.options = options;
out.verbose = verbose;

X = VBA_vec(X);
if isempty(options.mask)
    L = length(X);
else
    L = length(find(options.mask==1));
end

% estimate smoothness
if ~isempty(options.R)
    out.fwhm = RFT_smoothness(options.R);
else
    out.fwhm = RFT_smoothness(X);
end
OUTSTR{2} = ['Search volume: L=',num2str(L)];
OUTSTR{3} = ['Estimated smoothness: FWHM=',num2str(out.fwhm,'%3.1f')];
OUTSTR{4} = ['Number of resels: R=',num2str(L./out.fwhm,'%3.1f')];
if verbose
    disp(OUTSTR{2})
    disp(OUTSTR{3})
    disp(OUTSTR{4})
end

% estimate critical height threshold (peak-level)
gridu = 0:1e-3:10;
pgrid = RFT_Pval(gridu',0,1,out.fwhm,L,options.type,options.dof);
d = abs(options.FWER-pgrid);
out.xc = gridu(find(d==min(d)));
OUTSTR{5} = ['Critical height threshold [peak-level]: X>',num2str(out.xc,'%3.2f'),' (test size: FWER=',num2str(round(options.FWER*100)),'%)'];
if verbose
    disp(OUTSTR{5})
end

% peak- level inference
peaks.ind = RFT_localmax(X);
peaks.val = X(peaks.ind);
peaks.prft = RFT_Pval(peaks.val,0,1,out.fwhm,L,options.type,options.dof);
switch options.type
    case 'norm'
        peaks.punc = 1-VBA_spm_Ncdf(peaks.val,0,1);
    case 't'
        peaks.punc = 1-VBA_spm_Tcdf(peaks.val,options.dof);
    case 'F'
        peaks.punc = 1-VBA_spm_Fcdf(peaks.val,options.dof(1),options.dof(2));
end

% cluster-level inference
[clusters.ind,clusters.imax] = RFT_clusters(X,options.u,0);
nc = length(clusters.ind);
clusters.k = [];
clusters.prft = [];
for i=1:nc
    clusters.k(i) = length(clusters.ind{i});
    clusters.prft(i) = RFT_Pval(options.u,clusters.k(i),1,out.fwhm,L,options.type,options.dof);
end

% set-level inference
if nc>0
    bigC = find(clusters.k>=options.k);
else
    bigC = [];
end
set.c = length(bigC);
set.prft = RFT_Pval(options.u,options.k,set.c,out.fwhm,L,options.type,options.dof);

% E[number of voxel per cluster] and E[number of clusters]
out.Em = RFT_expectedTopo(options.u,L,out.fwhm,1,options.type,options.dof);
switch options.type
    case 'norm'
        P0 = 1-VBA_spm_Ncdf(options.u,0,1);
    case 't'
        P0 = 1-VBA_spm_Tcdf(options.u,options.dof);
    case 'F'
        P0 = 1-VBA_spm_Fcdf(options.u,options.dof(1),options.dof(2));
end
out.En = L.*P0./out.Em;

OUTSTR{6} = ['Expected voxels per cluster [cluster-level]: E[k|H0]=',num2str(out.En,'%3.1f')];
OUTSTR{7} = ['Expected number of clusters [set-level]: E[c|H0]=',num2str(out.Em,'%3.1f')];
if verbose
    disp(OUTSTR{6})
    disp(OUTSTR{7})
end

OUTSTR{8} = ['Number of local peaks =',num2str(length(peaks.prft)),];
OUTSTR{9} = ['Number of upcrossing clusters =',num2str(nc),' (X>',num2str(options.u,'%3.2f'),')'];
OUTSTR{10} = ['RFT [set-level]: p=',num2str(set.prft,'%3.3f'),' (c=',num2str(set.c),')'];
if verbose
    disp(OUTSTR{8})
    disp(OUTSTR{9})
    disp(OUTSTR{10})
end



% wrap-up results
out.peaks = peaks;
out.clusters = clusters;
out.set = set;
out.OUTSTR = OUTSTR;


% order peak-pval of clusters' maxima
inc = ismember(peaks.ind,clusters.imax);
[ps,is] = sort(peaks.prft(inc==1),'ascend');
STR.loc = cell(0);
STR.unc = cell(0);
STR.peak = cell(0);
STR.cluster = cell(0);
for i=1:nc % for all clusters
    % find local maxima that belong to each cluster
    ic = is(i);
    ip = find(ismember(peaks.ind,clusters.ind{ic})==1);
    % sort peak-pvalues in each cluster
    [pval,iop] = sort(peaks.prft(ip),'ascend');
    % aggregate info Re: local maxima per (ordered) cluster
    nt = length(STR.loc);
    for j=1:length(ip)
        STR.loc{nt+j} = num2str(peaks.ind(ip(iop(j))));
        STR.unc{nt+j} = num2str(peaks.punc(ip(iop(j))),'%3.3f');
        STR.peak{nt+j} = [num2str(peaks.prft(ip(iop(j))),'%3.3f'),'   (',num2str(peaks.val(ip(iop(j))),'%3.2f'),')'];
        if j==1
            STR.cluster{nt+j} = [num2str(clusters.prft(ic),'%3.3f'),'   (',num2str(clusters.k(ic)),')'];
        end
    end
end
STR.set = [num2str(set.prft,'%3.3f'),'   (',num2str(set.c),')'];
out.STR = STR;

% display results
if verbose
    [out] = RFT_ReDisplay(X,out);
end

