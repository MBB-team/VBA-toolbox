function [pv,stat,df,all] = GLM_contrast(X,y,c,type,verbose,Xnames,Ynames)
% computes classical p-values for any contrast applied onto GLM effects
% function [pv,stat,df] = GLM_contrast(X,y,c,type)
% IN:
%   - X: nXk design matrix
%   - y: nXp data matrix
%   - c: kXm contrast matrix
%   - type: flag for t- or F- test. Can be set to 't' (default) or 'F'
%   - verbose: flag for displaying results
%   - names: kx1 cell array of independent variables names
% OUT:
%   - pv: pX1 vector of p-values
%   - stat: pX1 vector of t- or F- statistics
%   - df: effective degees of freedom of the t- or F- test
%   - all: structure array with fields:
%       .R2: pX1 vector of coefficients of determination
%       .R2_a: pX1 vector of adjusted coefficients of determination (only
%       for F-test).
%       .b: kXp matrix of parameter OLS-estimates
%       .iC: kXk unscaled parameter covariance matrix
%       .vhat: kX1 vector of residal variances
%   NB: the covariance matrix Q of the i^th set of parameter is defined by:
%   Q = all.vhat(i).*all.iC.
%   If verbose=1, then 'all' also contains summary statistics of F-tests
%   wrt each parameter, through the following fields:
%       .pv: kXp matrix of p-values
%       .stat: kXp matrix of F- statistics
%   	.df: kX2Xp matrix of degees of freedom


[n,p] = size(y);
k = size(X,2);
try;type;catch;type='t';end
try;verbose;catch;verbose=0;end
try;Xnames{k};catch;Xnames=[];end
try;Ynames{p};catch;Ynames=[];end
C = X'*X;
iC = pinv(C);
b = iC*X'*y;
P = X*iC*X';
yhat = P*y;
R = eye(n) - P;
trR = trace(R);
stat = zeros(p,1);
pv = zeros(p,1);
vhat = zeros(p,1);
R2 = zeros(p,1);
switch type
    
    case 't'
        
        df = trR.^2./trace(R*R);
        for i=1:p
            vhat(i) = y(:,i)'*R*y(:,i)./trR;
            V = vhat(i).*c'*iC*c;
            stat(i) = c'*b(:,i)./sqrt(V);
            pv(i) = 1 - spm_Tcdf(stat(i),df);
            SS_tot = sum((y(:,i)-mean(y(:,i))).^2);
            SS_err = sum((y(:,i)-yhat(:,i)).^2);
            R2(i) = 1-(SS_err/SS_tot);
        end
        
    case 'F'
        
        ic = pinv(c'*c)*c';
        c0 = eye(size(c,1)) - c*ic;
        X0 = X*c0;
        R0 = eye(n) - X0*pinv(X0'*X0)*X0';
        y_a = R0*y;
        yhat_a = R0*yhat;
        R2_a = zeros(p,1);
        M = R0 - R;
        df = [trace(M).^2./trace(M*M),trR.^2./trace(R*R)];
        for i=1:p
            vhat(i) = y(:,i)'*R*y(:,i)./trR;
            stat(i) = ((b(:,i)'*X'*M*X*b(:,i))./(y(:,i)'*R*y(:,i))).*(trR./trace(R0-R));
            pv(i) = 1 - spm_Fcdf(stat(i),df(1),df(2));
            SS_tot = sum((y(:,i)-mean(y(:,i))).^2);
            SS_err = sum((y(:,i)-yhat(:,i)).^2);
            R2(i) = 1-(SS_err/SS_tot);
            SS_tot_a = sum((y_a(:,i)-mean(y_a(:,i))).^2);
            SS_err_a = sum((y_a(:,i)-yhat_a(:,i)).^2);
            R2_a(i) = 1-(SS_err_a/SS_tot_a);
        end
        
    otherwise
        
        disp('Error: this function only supports t- and F- tests!')
        pv = [];
        stat = [];
        df = [];
        all = [];
        return;
        
end

all.R2 = R2;
if isequal(type,'F')
    all.R2_a = R2_a;
end
all.b = b;
all.iC = iC;
all.vhat = vhat;

if ~verbose
    return;
end

% first run F-test through all regressors
all.pv = zeros(k,p);
all.stat = zeros(k,p);
all.df = zeros(k,2,p);
for i=1:p
    for j=1:k
        c = zeros(k,1);
        c(j) = 1;
        [all.pv(j,i),all.stat(j,i),all.df(j,:,i)] = GLM_contrast(X,y(:,i),c,'F',0);
    end
end

% axes for predicted vs observed data
handles.hf = figure('color',[1 1 1]);
handles.ha = subplot(2,2,1,...
    'parent',handles.hf,...
    'nextplot','add',...
    'visible','on');

% axes for parameter estimates
handles.ha(2) = subplot(2,2,2,...
    'parent',handles.hf,...
    'nextplot','add',...
    'visible','on');


% parameters' correlation matrix
handles.ha(3) = subplot(2,2,3,...
    'parent',handles.hf,...
    'nextplot','add');
imagesc(cov2corr(iC),'parent',handles.ha(3))
axis(handles.ha(3),'square')
axis(handles.ha(3),'equal')
axis(handles.ha(3),'tight')
colorbar('peer',handles.ha(3))
colormap(cool)
title(handles.ha(3),'parameters'' correlation matrix')
set(handles.ha(3),...
    'clim',[-1,1],...
    'xdir','normal',...
    'ydir','reverse',...
    'xtick',[1:k],...
    'ytick',[1:k])

% display test results (p-value and F-test)
ha = subplot(2,2,4,'parent',handles.hf,'visible','off');
pos = get(ha,'position');
un = get(ha,'units');
delete(ha)
dy = pos(4)/4;
pos1 = pos + [0,-3*dy+pos(4),0,dy-pos(4)];
handles.ht(1) = uicontrol(...
    'parent',handles.hf,...
    'style','text',...
    'units',un,...
    'position',pos1,...
    'string',[],...
    'fontsize',10,...
    'backgroundcolor',[1 1 1],...
    'HorizontalAlignment','left');

pos1 = pos + [0,-4*dy+pos(4),0,dy-pos(4)];
if length(df)==1
    str = ['dof=',num2str(df)];
else
    str = ['dof=[',num2str(df(1)),',',num2str(df(2)),']'];
end
handles.ht(2) = uicontrol(...
    'parent',handles.hf,...
    'style','text',...
    'units',un,...
    'position',pos1,...
    'string',str,...
    'fontsize',10,...
    'backgroundcolor',[1 1 1],...
    'HorizontalAlignment','left');

% data selector
str = cell(p,1);
for i=1:p
    str{i} = ['data #',num2str(i)];
    if ~isempty(Ynames)
        str{i} = [str{i},' (',Ynames{i},')'];
    end
end
ud.type = type;
ud.all = all;
ud.y = y;
ud.yhat = yhat;
ud.pv = pv;
ud.stat = stat;
ud.Xnames = Xnames;
ud.Ynames = Ynames;
ud.type = type;
ud.handles = handles;
pos1 = pos + [0,-2*dy+pos(4),0,dy-pos(4)];
handles.ht(2) = uicontrol(...
    'parent',handles.hf,...
    'style','popupmenu',...
    'units',un,...
    'position',pos1,...
    'string',str,...
    'fontsize',10,...
    'backgroundcolor',0.8*[1 1 1],...
    'HorizontalAlignment','left',...
    'userdata',ud,...
    'callback',@myData);
feval(@myData,handles.ht(2),[])



try,getSubplots;end



function colors = getColors(n)
hf = figure('visible','off');
ha = axes('parent',hf);
colors = get(ha,'colororder');
colors = repmat(colors,ceil(n/size(colors,1)),1);
colors = colors(1:n,:);
delete(hf)


function myData(e1,e2)
ud = get(e1,'userdata');
ind = get(e1,'value');

% plot data
cla(ud.handles.ha(1))
hp = plot(ud.handles.ha(1),ud.y(:,ind),ud.yhat(:,ind),'k.');
mi = min([ud.y(:,ind);ud.yhat(:,ind)]);
ma = max([ud.y(:,ind);ud.yhat(:,ind)]);
plot(ud.handles.ha(1),[mi,ma],[mi,ma],'r')
xx = 0.1*(ma-mi)+mi;
str = ['R^2=',sprintf('%2.3f',ud.all.R2(ind))];
if isequal(ud.type,'F')
    str = [str, ' [adj.R^2=',sprintf('%2.3f',ud.all.R2_a(ind)),']'];
end
text(xx,xx,str,'parent',ud.handles.ha(1),'color','r')
axis(ud.handles.ha(1),'tight')
xlabel(ud.handles.ha(1),'observed data')
ylabel(ud.handles.ha(1),'predicted data')
grid(ud.handles.ha(1),'on')
title(ud.handles.ha(1),'data fit')

% parameter estimates
cla(ud.handles.ha(2))
Vb = diag(ud.all.vhat(ind)*ud.all.iC);
k = size(ud.all.b,1);
for j=1:k
    hp = bar(ud.handles.ha(2),j,ud.all.b(j,ind),'facecolor',0.8*[1 1 1],'BarWidth',0.5);
    hcmenu = uicontextmenu;
    str = ['data #',num2str(ind)];
    if ~isempty(ud.Ynames)
        str = [str,' (',ud.Ynames{ind},')'];
    end
    str = [str,', variable #',num2str(j)];
    if ~isempty(ud.Xnames)
        str = [str,' (',ud.Xnames{j},')'];
    end
    uimenu(hcmenu, 'Label',str);
    uimenu(hcmenu, 'Label',['p=',num2str(ud.all.pv(j,ind),'%3.3f')]);
    uimenu(hcmenu, 'Label',['F=',num2str(ud.all.stat(j,ind),'%3.3f')]);
    uimenu(hcmenu, 'Label',['dof=[',num2str(ud.all.df(j,1,ind)),',',num2str(ud.all.df(j,2,ind)),']']);
    set(get(hp,'children'),'uicontextmenu',hcmenu);
    set(hp,'uicontextmenu',hcmenu);
    hp = errorbar(ud.handles.ha(2),j,ud.all.b(j,ind),1.96*sqrt(Vb(j)),'r.');
    set(hp,'uicontextmenu',hcmenu);
end
set(ud.handles.ha(2),'xtick',[1:1:k],'xlim',[0.5,k+0.5],'ygrid','on')
xlabel(ud.handles.ha(2),'independent variables')
title(ud.handles.ha(2),'parameter estimates')


% display test results
strp = ['p=',num2str(ud.pv(ind),'%3.3f'),' (',ud.type,'=',num2str(ud.stat(ind),'%3.3f'),')'];
set(ud.handles.ht(1),'string',strp);



function F = spm_Tcdf(x,v)
% Cumulative Distribution Function (CDF) of Students t distribution
% Copyright (C) 1992-2011 Wellcome Trust Centre for Neuroimaging
if nargin<2, error('Insufficient arguments'), end
ad = [ndims(x);ndims(v)];
rd = max(ad);
as = [[size(x),ones(1,rd-ad(1))];...
      [size(v),ones(1,rd-ad(2))]];
rs = max(as);
xa = prod(as,2)>1;
if all(xa) && any(diff(as(xa,:)))
    error('non-scalar args must match in size');
end
%-Initialise result to zeros
F = zeros(rs);
%-Only defined for strictly positive v. Return NaN if undefined.
md = ( ones(size(x))  &  v>0 );
if any(~md(:))
    F(~md) = NaN;
    warning('Returning NaN for out of range arguments');
end
%-Special case: f is 0.5 when x=0 (where betainc involves log of zero)
F( md  &  x==0 ) = 0.5;
%-Special case: Standard Cauchy distribution when v=1
ml = ( md  &  v==1 ); if xa(1), mlx=ml; else mlx=1; end
F(ml) = 0.5 + atan(x(mlx))/pi;
%-Compute where defined & not special cases
Q  = find( md  &  x~=0  &  v~=1 );
if isempty(Q), return, end
if xa(1), Qx=Q; else Qx=1; end
if xa(2), Qv=Q; else Qv=1; end
%-Compute
xQxPos = x(Qx)>0;
F(Q) = xQxPos -(xQxPos*2-1).*0.5.*betainc(v(Qv)./(v(Qv)+x(Qx).^2),v(Qv)/2,1/2);


function F = spm_Fcdf(x,v,w)
% Cumulative Distribution Function (CDF) of F (Fisher-Snedecor) distribution
% Copyright (C) 1992-2011 Wellcome Trust Centre for Neuroimaging
if nargin<2, error('Insufficient arguments'), end
%-Unpack degrees of freedom v & w from single df parameter (v)
if nargin<3
    vs = size(v);
    if prod(vs)==2
        %-DF is a 2-vector
        w = v(2); v = v(1);
    elseif vs(end)==2
        %-DF has last dimension 2 - unpack v & w
        nv = prod(vs);
        w  = reshape(v(nv/2+1:nv),vs(1:end-1));
        v  = reshape(v(1:nv/2)   ,vs(1:end-1));
    else
        error('Can''t unpack both df components from single argument')
    end
end
%-Check argument sizes
ad = [ndims(x);ndims(v);ndims(w)];
rd = max(ad);
as = [[size(x),ones(1,rd-ad(1))];...
      [size(v),ones(1,rd-ad(2))];...
      [size(w),ones(1,rd-ad(3))]];
rs = max(as);
xa = prod(as,2)>1;
if sum(xa)>1 && any(any(diff(as(xa,:)),1))
    error('non-scalar args must match in size'), end
%-Initialise result to zeros
F = zeros(rs);
%-Only defined for strictly positive v & w. Return NaN if undefined.
md = ( ones(size(x))  &  v>0  &  w>0 );
if any(~md(:))
    F(~md) = NaN;
    warning('Returning NaN for out of range arguments');
end
%-Non-zero where defined and x>0
Q  = find( md  &  x>0 );
if isempty(Q), return, end
if xa(1), Qx=Q; else Qx=1; end
if xa(2), Qv=Q; else Qv=1; end
if xa(3), Qw=Q; else Qw=1; end
%-Compute
F(Q) = 1 - betainc(w(Qw)./(w(Qw) + v(Qv).*x(Qx)),w(Qw)/2,v(Qv)/2);

