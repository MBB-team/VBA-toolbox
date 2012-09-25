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
try; names;catch;names = cell(p,1);end
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
all.pv = zeros(k,1);
all.stat = zeros(k,1);
all.df = zeros(k,2);
for j=1:k
    c = zeros(k,1);
    c(j) = 1;
    [all.pv(j),all.stat(j),all.df(j,:)] = GLM_contrast(X,y,c,'F',0);
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
for i=1:p
    Vb = diag(vhat(i)*iC);
    for j=1:k
        hp = bar(ha,j,b(j,i));
        hcmenu = uicontextmenu;
        uimenu(hcmenu, 'Label',['#',num2str(j),': ',names{j}]);
        uimenu(hcmenu, 'Label',['p=',num2str(all.pv(j))]);
        uimenu(hcmenu, 'Label',['F=',num2str(all.stat(j))]);
        uimenu(hcmenu, 'Label',['dof=[',num2str(all.df(j,1)),',',num2str(all.df(j,2)),']']);
        set(get(hp,'children'),'uicontextmenu',hcmenu);
        set(hp,'uicontextmenu',hcmenu);
        hp = errorbar(ha,j,b(j,i),1.96*sqrt(Vb(j)),'.','color',col(i,:));
        set(hp,'uicontextmenu',hcmenu);
    end
end
set(ha,'xtick',[1:1:k],'xlim',[0.5,k+0.5],'ygrid','on')
xlabel('independent variables')
title('parameter estimates')

% parameters' correlation matrix
ha = subplot(2,2,3,'parent',hf,'nextplot','add');
imagesc(cov2corr(iC),'parent',ha)
axis(ha,'square')
axis(ha,'equal')
axis(ha,'tight')
colorbar('peer',ha)
title(ha,'parameters'' correlation matrix')
set(ha,'clim',[-1,1])

ha = subplot(2,2,4,'parent',hf,'visible','off');
pos = get(ha,'position');
un = get(ha,'units');
delete(ha)
str{1} = ['p=',num2str(pv)];
switch type
    case 't'
        str{2} = ['t=',num2str(stat)];
        str{3} = ['dof=',num2str(df)];
    case 'F'
        str{2} = ['F=',num2str(stat)];
        str{3} = ['dof=[',num2str(df(1)),',',num2str(df(2)),']'];
end
hu = uicontrol(...
    'style','text',...
    'units',un,...
    'position',pos,...
    'string',str,...
    'fontsize',10,...
    'backgroundcolor',[1 1 1]);

try,getSubplots;end



