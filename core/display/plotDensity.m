function [h] = plotDensity(f_fname,g_fname,u,n_t,options,dim,pX,gX,pY,gY)

% plot predictive density 
% function [hp] =
% plotDensity(f_fname,g_fname,u,n_t,options,dim,pX,gX,pY,gY)
% IN:
%   - f_fname: name/handle of the evolution function
%   - g_fname: name/handle of the observation function
%   - u: input to the system
%   - n_t: number of time samples of the predictive density
%   - options: the optional structure. The prior pdfs on
%   evolution/observation parameters/hyperparameters are extracted from
%   this structure.
%   - dim: the dimension of the sDCM generative model
%   - pX: the n_tXnpXdim.n 3D array of MCMC empirical histograms of hidden
%   states
%   - gX: the npXdim.n 2D array giving the grid used for forming the MCMC
%   empirical histograms on each dimension of the hidden states
%   - pX/gY: [id, but for observed data]
% See also: VBA_MCMC_predictiveDensity.m

% fill in option structure
try
    priors0 = options.priors;
    if isinf(options.priors.a_alpha)
        options.priors.a_alpha = 1;
        options.priors.b_alpha = 1;
    end
end
options.verbose = 0;
dim.n_t = n_t;
[options] = VBA_check([],u,f_fname,g_fname,dim,options);
try
    options.priors.a_alpha = priors0.a_alpha;
    options.priors.b_alpha = priors0.b_alpha;
end

% get deterministic trajectory
[y,x] = VBA_simulate(...
    n_t,...
    f_fname,...
    g_fname,...
    options.priors.muTheta,...
    options.priors.muPhi,...
    u,...
    Inf,...
    Inf,...
    options,...
    options.priors.muX0);
    
hfp = figure(...
        'color',[1 1 1],...
        'name','densities',...
        'tag','densities',...
        'menu','none',...
        'Renderer','OpenGL');
labels = {'summary','observables','states'};
callbacks = {@mySummary,@myObs,@myStates};
[h] = VBA_spm_uitab(hfp,labels,callbacks,'XY',1,1,0.05);

h.hfp = hfp;
in = struct(...
    'h',h,...
    'pX',pX,...
    'gX',gX,...
    'pY',pY,...
    'gY',gY,...
    'y',y,...
    'x',x,...
    'options',options,...
    'dim',dim);
set(hfp,'userdata',in)
set(h.htab,'backgroundcolor',[1 1 1])
set(h.hh,'backgroundcolor',[1 1 1])
set(h.hp,'HighlightColor',0.8*[1 1 1])
set(h.hp,'backgroundcolor',[1 1 1])
mySummary(hfp)


function [h] = mySummary(hf)
try 
    in = get(hf,'userdata');
catch
    in = get(gcf,'userdata');
end
try; delete(get(in.h.hp,'children')); end
options = in.options;
dim = in.dim;
str{1} = sprintf(['Dimensions of the model:','\n ',...
    '    - data: p=',num2str(dim.p),'\n ',...
    '    - time samples: t=',num2str(dim.n_t),'\n ',...
    '    - hidden states: n=',num2str(dim.n),'\n ',...
    '    - evolution parameters: n_theta=',num2str(dim.n_theta),'\n ',...
    '    - observation parameters: n_phi=',num2str(dim.n_phi),'\n ']);
if numel(options.sources) > 1
    tmp = ' (multisources)';
else
    switch options.sources.type
        case 0
            tmp = ' (gaussian data)';
        case 1
            tmp = ' (binomial data)';
        case 2
            tmp = ' (multinomial data)';
    end
end
if dim.n >= 1
    if isinf(options.priors.a_alpha) ...
            && isequal(options.priors.b_alpha,0)
        str{2} = 'This was a deterministic dynamical system';
    else
        str{2} = 'This was a stochastic dynamical system';
    end
    if isa(options.g_fname,'function_handle')
        gfn = func2str(options.g_fname);
    else
        gfn = options.g_fname;
    end
    if isequal(gfn,'g_embed')
        gfn0 = options.inG.g_fname;
        if isa(gfn0,'function_handle')
            gfn0 = func2str(gfn0);
        end
        gfn = [gfn,' (',gfn0,')'];
        str{2} = [str{2},' (with delay embedding)'];
    end
    if isa(options.f_fname,'function_handle')
        ffn = func2str(options.f_fname);
    else
        ffn = options.f_fname;
    end
    if isequal(ffn,'f_embed')
        ffn0 = options.inF.f_fname;
        if isa(ffn0,'function_handle')
            ffn0 = func2str(ffn0);
        end
        ffn = [ffn,' (',ffn0,')'];
    end
    str{3} = sprintf([...
        '    - observation function: ',gfn,tmp,'\n',...
        '    - evolution function: ',ffn,'\n ']);
else
    str{2} = 'The model was static (no hidden states)';
    if isa(options.g_fname,'function_handle')
        gfn = func2str(options.g_fname);
    else
        gfn = options.g_fname;
    end
    str{5} = sprintf(['    - observation function: ',gfn,tmp,'\n ']);
end
uicontrol(...
    'parent',in.h.hp,...
    'style','text',...
    'tag','plot',...
    'units','normalized',...
    'position',[0.1,0.1,0.8,0.8],...
    'backgroundcolor',get(in.h.hp,'BackgroundColor'),...
    'HorizontalAlignment','left',...
    'fontsize',11,...
    'string',str);


function [h] = myObs(ho,evt)
in = get(get(gco,'parent'),'userdata');
try; delete(get(in.h.hp,'children')); end
n = size(in.y,1);
for i=1:n
   labels{i} = ['dim ',num2str(i)];
   callbacks{i} = @plotPY;
end
[h] = VBA_spm_uitab(in.h.hp,labels,callbacks,'plot',1,1,0.05);
set(h.htab,'backgroundcolor',[1 1 1])
set(h.hh,'backgroundcolor',[1 1 1])
set(h.hp,'HighlightColor',0.8*[1 1 1])
set(h.hp,'backgroundcolor',[1 1 1])
plotPY(h.htab(1))


function [h] = myStates(ho,evt)
in = get(get(gco,'parent'),'userdata');
try; delete(get(in.h.hp,'children')); end
n = size(in.x,1);
for i=1:n
   labels{i} = ['dim ',num2str(i)];
   callbacks{i} = @plotPX;
end
[h] = VBA_spm_uitab(in.h.hp,labels,callbacks,'plot',1,1,0.05);
set(h.htab,'backgroundcolor',[1 1 1])
set(h.hh,'backgroundcolor',[1 1 1])
set(h.hp,'HighlightColor',0.8*[1 1 1])
set(h.hp,'backgroundcolor',[1 1 1])
plotPX(h.htab(1))


function [h] = plotPX(ho)
try 
    in0 = get(ho,'userdata');
catch
    in0 = get(gco,'userdata');
end
in = get(gcf,'userdata');
try; delete(get(in0.handles.hp,'children')); end
ha = axes('parent',in0.handles.hp,'tag','graph');
ud = struct('ind',in0.ind,'ha',ha,'XY','x');
uicontrol(...
    'parent',in0.handles.hp,...
    'style','edit',...
    'tag','graph',...
    'units','normalized',...
    'backgroundcolor',[1 1 1],...
    'HorizontalAlignment','left',...
    'fontsize',11,...
    'userdata',ud,...
    'tooltipstring','variance of the smoothing kernel',...
    'callback',@resmooth,...
    'string','1');
plotGraph3D(in.pX(:,:,in0.ind),in.gX(:,in0.ind),in.x(in0.ind,:),ha,1);


function [h] = plotPY(ho)
try 
    in0 = get(ho,'userdata');
catch
    in0 = get(gco,'userdata');
end
in = get(gcf,'userdata');
try; delete(get(in0.handles.hp,'children')); end
ha = axes('parent',in0.handles.hp,'tag','graph');
ud = struct('ind',in0.ind,'ha',ha,'XY','y');
uicontrol(...
    'parent',in0.handles.hp,...
    'style','edit',...
    'tag','graph',...
    'units','normalized',...
    'backgroundcolor',[1 1 1],...
    'HorizontalAlignment','left',...
    'fontsize',11,...
    'userdata',ud,...
    'tooltipstring','variance of the smoothing kernel',...
    'callback',@resmooth,...
    'string','1');
plotGraph3D(in.pY(:,:,in0.ind),in.gY(:,in0.ind),in.y(in0.ind,:),ha,1);


function resmooth(ho,evt)
ud = get(gco,'userdata');
va = abs(str2double(get(gco,'string')));
in = get(gcf,'userdata');
switch ud.XY
    case 'x'
        plotGraph3D(in.pX(:,:,ud.ind),in.gX(:,ud.ind),in.x(ud.ind,:),ud.ha,va);
    case 'y'
        plotGraph3D(in.pY(:,:,ud.ind),in.gY(:,ud.ind),in.y(ud.ind,:),ud.ha,va);
end

