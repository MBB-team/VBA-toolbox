function [pv,stat,df,all] = GLM_contrast(X,y,c,type,verbose,names)
% computes classical p-values for any contrast applied onto GLM effects
% function [pv,stat,df] = GLM_contrast(X,y,c,type)
% IN:
%   - X: nXk designa matrix
%   - y: nXp data matrix
%   - c: kXm contrast matrix
%   - type: flag for t- or F- test. Can be set to 't' (default) or 'F'
%   - verbose: flag for displaying results
%   - names: kx1 cell array of independent variables names
% OUT:
%   - pv: pX1 vector of p-values
%   - stat: pX1 vector of t- or F- statistics
%   - df: effective degees of freedom of the t- or F- test
%   - all: structure array with fields .p, .stat and .df, which are summary
%   statistics of F-tests for each and every independent variable included
%   in the analysis.


[n,p] = size(y);
k = size(X,2);
try;type;catch;type='t';end
try;verbose;catch;verbose=0;end
try; names{k};catch;names = cell(k,1);end
C = X'*X;
iC = pinv(C);
b = iC*X'*y;
P = X*iC*X';
R = eye(n) - P;
trR = trace(R);
stat = zeros(p,1);
pv = zeros(p,1);
vhat = zeros(p,1);

switch type
    
    case 't'
        
        df = trR.^2./trace(R*R);
        for i=1:p
            vhat(i) = y(:,i)'*R*y(:,i)./trR;
            V = vhat(i).*c'*iC*c;
            stat(i) = c'*b(:,i)./sqrt(V);
            pv(i) = 1 - spm_Tcdf(stat(i),df);
        end
        
    case 'F'
        
        ic = pinv(c'*c)*c';
        c0 = eye(size(c,1)) - c*ic;
        X0 = X*c0;
        R0 = eye(n) - X0*pinv(X0'*X0)*X0';
        M = R0 - R;
        df = [trace(M).^2./trace(M*M),trR.^2./trace(R*R)];
        for i=1:p
            vhat(i) = y(:,i)'*R*y(:,i)./trR;
            stat(i) = ((b(:,i)'*X'*M*X*b(:,i))./(y(:,i)'*R*y(:,i))).*(trR./trace(R0-R));
            pv(i) = 1 - spm_Fcdf(stat(i),df(1),df(2));
        end
        
    otherwise
        
        disp('Error: this function only supports t- and F- tests!')
        out = [];
        
end


if ~verbose
    all = [];
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


col = getColors(p);

% data fit
yhat = P*y;
hf = figure('color',[1 1 1]);
ha = subplot(2,2,1,'parent',hf,'nextplot','add');
for i=1:p
    hp = plot(ha,y(:,i),yhat(:,i),'.','color',col(i,:));
end
mi = min([y(:);yhat(:)]);
ma = max([y(:);yhat(:)]);
plot([mi,ma],[mi,ma],'r')
axis(ha,'tight')
xlabel(ha,'observed data')
ylabel(ha,'predicted data')
grid(ha,'on')
title(ha,'data fit')

% parameter estimates
ha = subplot(2,2,2,'parent',hf,'nextplot','add');
offset = [1:p]/(2*p);
offset = offset - mean(offset);
width = diff(offset);
if isempty(width)
    width = 0.5;
else
    width = width(1);
end
for i=1:p
    Vb = diag(vhat(i)*iC);
    for j=1:k
        hp = bar(ha,j+offset(i),b(j,i),'facecolor',col(i,:),'BarWidth',width);
        hcmenu = uicontextmenu;
        uimenu(hcmenu, 'Label',['data #',num2str(i),', variable #',num2str(j),' (',names{j},')']);
        uimenu(hcmenu, 'Label',['p=',num2str(all.pv(j,i),'%3.3f')]);
        uimenu(hcmenu, 'Label',['F=',num2str(all.stat(j,i),'%3.3f')]);
        uimenu(hcmenu, 'Label',['dof=[',num2str(all.df(j,1,i)),',',num2str(all.df(j,2,i)),']']);
        set(get(hp,'children'),'uicontextmenu',hcmenu);
        set(hp,'uicontextmenu',hcmenu);
        hp = errorbar(ha,j+offset(i),b(j,i),1.96*sqrt(Vb(j)),'.','color',col(i,:));
        set(hp,'uicontextmenu',hcmenu);
    end
end
set(ha,'xtick',[1:1:k],'xlim',[0.5+min(offset),k+0.5+max(offset)],'ygrid','on')
xlabel('independent variables')
title('parameter estimates')

% parameters' correlation matrix
ha = subplot(2,2,3,'parent',hf,'nextplot','add');
imagesc(cov2corr(iC),'parent',ha)
axis(ha,'square')
axis(ha,'equal')
axis(ha,'tight')
colorbar('peer',ha)
colormap(cool)
title(ha,'parameters'' correlation matrix')
set(ha,...
    'clim',[-1,1],...
    'xdir','normal',...
    'ydir','reverse',...
    'xtick',[1:k],...
    'ytick',[1:k])

% display test results
ha = subplot(2,2,4,'parent',hf,'visible','off');
pos = get(ha,'position');
un = get(ha,'units');
delete(ha)
for i=1:p
    dy = i*pos(4)/8;
    pos1 = pos + [0,-dy+pos(4),0,dy/i-pos(4)];
    strp = ['data #',num2str(i),': p=',num2str(pv(i),'%3.3f'),' (',type,'=',num2str(stat(i),'%3.3f'),')'];
    hj = uicontrol(...
        'style','text',...
        'units',un,...
        'position',pos1,...
        'string',strp,...
        'fontsize',10,...
        'backgroundcolor',[1 1 1],...
        'ForegroundColor',col(i,:),...
        'HorizontalAlignment','left');
end
if length(df)==1
    str = ['dof=',num2str(df)];
else
    str = ['dof=[',num2str(df(1)),',',num2str(df(2)),']'];
end
pos2 = pos1 + [0,-dy/i,0,0];
hu = uicontrol(...
    'style','text',...
    'units',un,...
    'position',pos2,...
    'string',str,...
    'fontsize',10,...
    'backgroundcolor',[1 1 1],...
    'HorizontalAlignment','left');



try,getSubplots;end



function colors = getColors(n)
hf = figure('visible','off');
ha = axes('parent',hf);
colors = get(ha,'colororder');
colors = repmat(colors,ceil(n/size(colors,1)),1);
colors = colors(1:n,:);
delete(hf)



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

