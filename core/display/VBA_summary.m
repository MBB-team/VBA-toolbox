function str = VBA_summary(out,newlines)
% writes a summary string from standard output of VBA model inversion
% function str = VBA_summary(out)
% IN:
%   - out: the 'out' structure of VBA inversion routine
%   - newlines: flag for inserting line separators (\n)
% OUT:
%   - str: a cell array of strings, which summarize the VBA inversion

try;newlines;catch,newlines=0;end

[LLH0] = VBA_LMEH0(out.y,out.options);
try F = out.F(end); catch, F = '?'; end
many = length(out.options.sources)>1; % more than one source?

str{1} = sprintf(['Date: ',datestr(out.date)]);
if ~out.options.OnLine
    s0 = ['VB converged in ',num2str(out.it),' iterations'];
else
    s0 = ['Online VB algorithm'];
end
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
if many
    datastr = [' (number of sources=',num2str(length(out.options.sources)),')'];
else
    datastr = [];
end
str{3} = sprintf(['Dimensions of the model:','\n ',...
    '    - data: p=',num2str(out.dim.p),datastr,'\n ',...
    '    - time samples: t=',num2str(out.dim.n_t),'\n ',...
    '    - hidden states: n=',num2str(out.dim.n),'\n ',...
    '    - evolution parameters: n_theta=',num2str(out.dim.n_theta),'\n ',...
    '    - observation parameters: n_phi=',num2str(out.dim.n_phi),'\n ',...
    '    - inputs: n_u=',num2str(out.dim.u)]);
if numel(out.options.sources)>1
    tmp = ' (multisource)';
else
    switch out.options.sources.type
        case 0
            tmp = ' (gaussian data)';
        case 1
            tmp = ' (binomial data)';
        case 2
            tmp = ' (multinomial data)';
    end
end
if out.options.UNL
    so = 'un-normalized likelihood';
else
    so = 'observation';
end
if out.dim.n >= 1
    if isinf(out.options.priors.a_alpha) && isequal(out.options.priors.b_alpha,0)
        str{4} = sprintf('This was a deterministic dynamical system');
    else
        str{4} = sprintf('This was a stochastic dynamical system');
    end
    if isa(out.options.g_fname,'function_handle')
        gfn = func2str(out.options.g_fname);
    else
        gfn = out.options.g_fname;
    end
    if isequal(gfn,'g_embed')
        gfn0 = out.options.inG.g_fname;
        if isa(gfn0,'function_handle')
            gfn0 = func2str(gfn0);
        end
        gfn = [gfn,' (',gfn0,')'];
        str{4} = [str{4},' (with delay embedding)'];
    end
    if isa(out.options.f_fname,'function_handle')
        ffn = func2str(out.options.f_fname);
    else
        ffn = out.options.f_fname;
    end
    if isequal(ffn,'f_embed')
        ffn0 = out.options.inF.f_fname;
        if isa(ffn0,'function_handle')
            ffn0 = func2str(ffn0);
        end
        ffn = [ffn,' (',ffn0,')'];
    end
    str{4} = sprintf([str{4},'\n ',...
        '    - ',so,' function: ',gfn,tmp,'\n ',...
        '    - evolution function: ',ffn]);
else
    str{4} = ['The model was static (no hidden states)','\n '];
    if isa(out.options.g_fname,'function_handle')
        gfn = func2str(out.options.g_fname);
    else
        gfn = out.options.g_fname;
    end
    str{4} = sprintf([str{4},'    - ',so,' function: ',gfn,tmp]);
end
str{5} = sprintf(['Bayesian log model evidences:','\n',...
    '     - full model: log p(y|m) > ',num2str(F,'%4.3e'),'\n',...
    '     - null hypothesis: log p(y|H0) = ',num2str(LLH0,'%4.3e')]);
if ~out.options.OnLine && out.dim.n >= 1 && ~isinf(out.options.priors.a_alpha) && ~isequal(out.options.priors.b_alpha,0)
    Fd = out.options.init.out.F;
    str{5} = sprintf([str{5},'\n ',...
        '    - deterministic variant: log p(y|m,eta=0) > ',num2str(Fd,'%4.3e') ]);
end
% str{5} = [str{5}, '\n'];

gsi = find([out.options.sources.type]==0);
bsi = find([out.options.sources.type]~=0);
many = length(out.options.sources)>1;
R2str = ['     - determin. coeff. (R2):  '];
LLstr = ['     - log-likelihood: '];
AICstr = ['     - AIC: '];
BICstr = ['     - BIC: '];
if ~isempty(gsi)
    R2str = [R2str,catnum2str(out.fit.R2,gsi,many,'%')];
    CAstr = [];
    LLstr = [LLstr,catnum2str(out.fit.LL,gsi,many)];
    AICstr = [AICstr,catnum2str(out.fit.AIC,gsi,many)];
    BICstr = [BICstr,catnum2str(out.fit.BIC,gsi,many)];
    separator = ', ';
else
    separator = [];
end
if  ~isempty(bsi)
    R2str = [R2str,separator,catnum2str(out.fit.R2,bsi,many,'%')];
    CAstr = ['     - balanced classif. acc.:  '];
    CAstr = [CAstr,catnum2str(out.fit.bacc,bsi,many,'%'),'\n'];
    LLstr = [LLstr,separator,catnum2str(out.fit.LL,bsi,many)];
    AICstr = [AICstr,separator,catnum2str(out.fit.AIC,bsi,many)];
    BICstr = [BICstr,separator,catnum2str(out.fit.BIC,bsi,many)];
end
R2str = [R2str,'\n'];
LLstr = [LLstr,'\n'];
AICstr = [AICstr,'\n'];
% BICstr = [BICstr,'\n'];
str{6} = ['Classical goodness-of-fit metrics:','\n',...
    R2str,...
    CAstr,...
    LLstr,...
    AICstr,...
    BICstr];

if newlines
    for i=1:length(str)
        str{i} = [str{i},'\n'];
    end
end


function str = catnum2str(x,ind,many,flag)
try,flag;catch,flag=[];end
str = [];
for i=1:length(ind)
    si=ind(i);
    if isequal(flag,'%')
        str = [str,', ',num2str(100*x(si),'%2.1f%%'),'%'];
    else
        str = [str,', ',num2str(x(si),'%4.3e')];
    end
    if many
        str = [str,' (source #',num2str(si),')'];
    end
end
str(1:2) = [];


