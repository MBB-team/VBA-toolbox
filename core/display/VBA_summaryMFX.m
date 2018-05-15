function str = VBA_summaryMFX(out)
% writes a summary string from output of VBA_MFX model inversion
% function str = VBA_summaryMFX(out)
% IN:
%   - out: the 'out' structure of VBA_MFX inversion routine
% OUT:
%   - str: a cell array of strings, which summarize the VBA inversion

str{1} = sprintf(['Date: ',datestr(out.date)]);
s0 = ['VB converged in ',num2str(out.it),' iterations'];
try
    if floor(out.dt./60) == 0
        timeString = [num2str(floor(out.dt)),' sec'];
    else
        timeString = [num2str(floor(out.dt./60)),' min'];
    end
    str{2} = sprintf([s0,' (took ~',timeString,')']);
catch
    str{2} = sprintf(s0);
end
str{3} = sprintf(['Dimensions of the MFX analysis:','\n ',...
    '    - subjects: ns=',num2str(out.options.dim.ns),'\n ',...
    '    - data: p=',num2str(out.options.dim.p),'\n ',...
    '    - hidden states: n=',num2str(out.options.dim.n),'\n ',...
    '    - evolution parameters: n_theta=',num2str(out.options.dim.n_theta),'\n ',...
    '    - observation parameters: n_phi=',num2str(out.options.dim.n_phi)]);
str{4} = sprintf('\n     Within-subject generative model:\n');
if isa(out.options.g_fname,'function_handle')
    gfn = func2str(out.options.g_fname);
else
    gfn = out.options.g_fname;
end
if out.options.dim.n >= 1
    if isa(out.options.f_fname,'function_handle')
        ffn = func2str(out.options.f_fname);
    else
        ffn = out.options.f_fname;
    end
    str{4} = sprintf([str{4},...
        '       - observation function: ',gfn,'\n',...
        '       - evolution function: ',ffn]);
else
    str{4} = sprintf([str{4},'     - observation function: ',gfn]);
end
mF = mean(out.within_fit.F(:));
sF = std(out.within_fit.F(:));
mF0 = mean(out.within_fit.LLH0(:));
sF0 = std(out.within_fit.LLH0(:));
str{5} = sprintf([...
    '\n     - Bayesian log model evidences: <log p(y|m)> = ',num2str(mF,'%4.3e'),' +/- ',num2str(sF,'%4.3e'),'\n',...
    '     - Bayesian log evidences under the null: <log p(y|H0)> = ',num2str(mF0,'%4.3e'),' +/- ',num2str(sF0,'%4.3e')]);

if any([out.options.sources.type]==0)
    R2str = 'coefficient of determination: <R2>';
    mR = mean(out.within_fit.R2);
    sR = std(out.within_fit.R2);
    str{6} = sprintf([...
        '\n     - ',R2str,' = ',num2str(mR,'%4.3f '),' +/- ',num2str(sR,'%4.3e ')]);
end

if any([out.options.sources.type]>0)
    R2str = 'balanced classification accuracy <Acc>';
    mR = mean(out.within_fit.R2);
    sR = std(out.within_fit.R2);
    str{7} = sprintf([...
        '\n     - ',R2str,' = ',num2str(mR,'%4.3f '),' +/- ',num2str(sR,'%4.3f ')]);
end


