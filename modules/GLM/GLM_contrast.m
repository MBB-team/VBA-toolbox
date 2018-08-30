function [pv,stat,df,all] = GLM_contrast(X,y,c,type,verbose,Xnames,Ynames,notest)
% computes classical p-values for any contrast applied onto GLM effects
% function [pv,stat,df] = GLM_contrast(X,y,c,type,verbose,Xnames,Ynames)
% In brief, this function uses classical tests of linear mixtures of
% effects, under the general linear model (GLM):
%   y = X*beta + e
% where y is the data, e is some (Gaussian) random noise, X is the design
% matrix and beta are unknown GLM parameters.
% Tests are specified in terms of a contrast matrix c, such that the
% p-value computes the probability of the test score under the null H0.
% More precisely:
%  - 't' tests are one-sided: H0={c'*beta>0}
%  - 'F' tests are two-sided: H0={c(:,1)'beta=0 AND ... c(:,m)'*beta=0}
% For example, testing for the ith parameter corresponds to an all-zero
% contrast vector, except on its ith entry.
% Note: To derive matricial F-contrasts for a main effect of an IV with n
% different values, see Contrast_MEbins.m!
% IN:
%   - X: nxk design matrix
%   - y: nxp data matrix
%   - c: kxm contrast matrix (default is eye(k) -> omnibus F-test)
%   - type: flag for t- or F- test. Can be set to 't' (default) or 'F'
%   - verbose: flag for displaying results (default is 0)
%   - Xnames: px1 cell array of independent variables names
%   - Ynames: kx1 cell array of dependent variables names
%   - notest: flag for testing all regressors in verbose mode
% OUT:
%   - pv: px1 vector of p-values
%   - stat: px1 vector of t- or F- statistics
%   - df: effective degees of freedom of the t- or F- test
%   - all: structure array with fields:
%       .ks: px1 vector of p-values for Kolmogorov-Smirnoff test
%       (residuals' normality)
%       .R2: px1 vector of coefficients of determination
%       .R2_a: px1 vector of amount of variance explained by the contrast
%       of interest
%       .b: kxp matrix of parameter OLS-estimates
%       .iC: kxk unscaled parameter covariance matrix
%       .vhat: px1 vector of residual variances
%   NB: the covariance matrix Q of the i^th set of parameter is defined by:
%   Q = all.vhat(i).*all.iC.
%   If verbose=1, then 'all' also contains summary statistics of F-tests
%   wrt each parameter, through the following fields:
%       .pv: kxp matrix of p-values
%       .stat: kxp matrix of F- statistics
%   	.df: kx2xp matrix of degees of freedom

% fill in default I/O
pv = [];
stat = [];
df = [];
all = [];
[n,p] = size(y);
[n0,k] = size(X);
[k0,m] = size(c);
try;c;catch;c=eye(k);type='F';end
try;type;catch;type='t';end
try;verbose;catch;verbose=0;end
try;Xnames{k};catch;Xnames=[];end
try;Ynames{p};catch;Ynames=[];end
try;notest;catch;notest=0;end

if verbose
    fprintf(1,'---')
    fprintf(1,'\n')
    fprintf(1,['Date: ',datestr(clock),'\n'])
    fprintf(1,'GLM contrast testing:\n')
end

% check basic numerical requirements
try
    if VBA_isWeird (y)
        disp('Error: data contains weird values!')
        return
    end
end
if ~isequal(n,n0)
    disp('Error: design matrix has to have as many rows as the data matrix!')
    return
end
if ~isequal(k,k0)
    disp('Error: contrasts have to have as many rows as there are columns in the design matrix!')
    return
end


% Estimate GLM parameters
if verbose
    fprintf(1,'Computing parameter covariance matrix...')
end
C = X'*X;
iC = pinv(C);
if verbose
    fprintf(1,' OK.')
    fprintf(1,'\n')
end
iX = iC*X'; % pseudo-inverse of design matrix X
b = iX*y;
if verbose
    fprintf(1,'Computing projection matrices...')
end
P = X*iX;
R = speye(n) - P;
if verbose
    fprintf(1,' OK.')
    fprintf(1,'\n')
end
yhat = X*b;
e = y-yhat;
trR = n - trace(P);

% perform significance testing
stat = zeros(p,1);
pv = zeros(p,1);
vhat = zeros(p,1);
R2 = zeros(p,1);
R2_a = zeros(p,1);

switch type
    
    case 't'
        
        if m~=1
            disp('Error: cannot have contrast matrices for ''t''-tests!')
            return
        end
        df = trR.^2./sum(sum(R.^2,1)); % Satterthwaite approx.
        if verbose
            fprintf(1,'Computing t-statistics...')
            if p>1
                fprintf(1,'%6.2f %%',0)
            end
        end
        vhat = sum((y-yhat).^2,1)./trR;
        V = vhat.*(c'*iC*c);
        stat = (c'*b)./sqrt(V);
        pv = 1 - VBA_spm_Tcdf(stat,df);
        for i=1:p
            SS_tot = sum((y(:,i)-mean(y(:,i))).^2);
            SS_err = sum(e(:,i).^2);
            R2(i) = 1-(SS_err/SS_tot);
            R2_a(i) = VBA_FtoR2(stat(i).^2,1,df);
            [tmp,ks(i)] = VBA_kstest(zscore(e));
            if verbose && p>1
                fprintf(1,repmat('\b',1,8))
                fprintf(1,'%6.2f %%',100*i/p)
            end
        end
        if verbose
            if p>1
                fprintf(1,repmat('\b',1,8))
            end
            fprintf(1,[' OK.'])
            fprintf(1,'\n')
        end
        
    case 'F'
        
        if verbose
            fprintf(1,'Computing contrast null-space projectors...')
        end
        ic = pinv(c'*c)*c';
        c0 = speye(size(c,1)) - c*ic;
        X0 = X*c0;
        R0 = speye(n) - X0*pinv(X0'*X0)*X0';
%         y_a = R0*y;
%         yhat_a = R0*yhat;
        M = R0 - R;
        trM = trace(M);
        df = [trM.^2./sum(sum(M.^2,1)),trR.^2./sum(sum(R.^2,1))];  % Satterthwaite approx.
        if verbose
            fprintf(1,' OK.')
            fprintf(1,'\n')
        end
        if verbose
            fprintf(1,'Computing F-statistics...')
            if p>1
                fprintf(1,'%6.2f %%',0)
            end
        end
        for i=1:p
            vhat(i) = sum(e(:,i).^2)./trR;
            stat(i) = ((yhat(:,i)'*M*yhat(:,i))./(y(:,i)'*R*y(:,i))).*(trR./trM);
            pv(i) = 1 - VBA_spm_Fcdf(stat(i),df(1),df(2));
            SS_tot = sum((y(:,i)-mean(y(:,i))).^2);
            SS_err = sum(e(:,i).^2);
            R2(i) = 1-(SS_err/SS_tot);
            R2_a(i) = VBA_FtoR2(stat(i),df(1),df(2));
%             SS_tot_a = sum((y_a(:,i)-mean(y_a(:,i))).^2);
%             SS_err_a = sum((y_a(:,i)-yhat_a(:,i)).^2);
%             R2_a(i) = 1-(SS_err_a/SS_tot_a);
            [tmp,ks(i)] = VBA_kstest(zscore(e));
            if verbose && p>1
                fprintf(1,repmat('\b',1,8))
                fprintf(1,'%6.2f %%',100*i/p)
            end
        end
        if verbose
            if p>1
                fprintf(1,repmat('\b',1,8))
            end
            fprintf(1,[' OK.'])
            fprintf(1,'\n')
        end
        
    otherwise
        
        disp('Error: this function only supports t- and F- tests!')
        return
        
end

% fill in output structure
all.R2 = R2;
all.R2_a = R2_a;
all.b = b;
all.iC = iC;
all.vhat = vhat;
all.ks = ks; % kolmogorov-smirnov test (normality of the residuals)
% all.tolerance = GLM_tolerance(X);

if ~verbose
    return;
end

% run F-test through all regressors
if ~notest
    all.pv = zeros(k,p);
    all.stat = zeros(k,p);
    all.df = zeros(k,2,p);
    if verbose
        fprintf(1,'Testing for each regressor significance...')
        if p*k>1
            fprintf(1,'%6.2f %%',0)
        end
    end
    for i=1:p
        for j=1:k
            cij = zeros(k,1);
            cij(j) = 1;
            [all.pv(j,i),all.stat(j,i),all.df(j,:,i)] = GLM_contrast(X,y(:,i),cij,'F',0);
            if verbose && p*k>1
                fprintf(1,repmat('\b',1,8))
                fprintf(1,'%6.2f %%',100*((i-1)*k+j)/(k*p))
            end
        end
    end
    if verbose
        if p*k>1
            fprintf(1,repmat('\b',1,8))
        end
        fprintf(1,' OK.')
        fprintf(1,'\n')
    end
end

% summarize results in matlab window
disp('COI significance and effect size:')
if length(df)==1
    strdf = ['dof=',num2str(df)];
else
    strdf = ['dof=[',num2str(df(1)),',',num2str(df(2)),']'];
end
stryn = cell(p,1);
for i=1:p
    stryn{i} = ['data #',num2str(i)];
    if ~isempty(Ynames)
        stryn{i} = [stryn{i},' (',Ynames{i},')'];
    end
    switch type
        case 't'
            disp([' - ',stryn{i},': p=',num2str(pv(i),3),', ',strdf])
        case 'F'
            disp([' - ',stryn{i},': p=',num2str(pv(i),3),', R2=',num2str(round(all.R2_a(i)*1e3)/10),'%, ',strdf])
    end
end

% create display figure
pos0 = get(0,'screenSize');
pos = [0.51*pos0(3),0.05*pos0(4),0.45*pos0(3),0.9*pos0(4)];
handles.hf = figure('color',[1 1 1],'position',pos,'menubar','none');

% axes for predicted vs observed data
handles.ha = subplot(3,2,1,'parent',handles.hf,'nextplot','add','visible','on');

% display design matrix
handles.ha(5) = subplot(3,2,6,'parent',handles.hf,'visible','on');
imagesc(X,'parent',handles.ha(5))
colorbar('peer',handles.ha(5))
set(handles.ha(5),'xtick',[1:1:k],'xlim',[0.5,k+0.5],'xgrid','on','ygrid','off')
xlabel(handles.ha(5),'independent variables')
ylabel(handles.ha(5),'dependent variables')
VBA_title(handles.ha(5),'design matrix')
pos = get(handles.ha(5),'position');

% axes for parameter estimates
handles.ha(2) = subplot(3,2,2,'parent',handles.hf,'nextplot','add','visible','on');
pos0 = get(handles.ha(2),'position');
set(handles.ha(2),'position',[pos(1),pos0(2),pos(3),pos0(4)])

% axes for predicted and observed data
handles.ha(6) = subplot(3,2,3,'parent',handles.hf,'nextplot','add','visible','on');

% parameters' correlation matrix
handles.ha(3) = subplot(3,2,5,'parent',handles.hf,'nextplot','add');
imagesc(VBA_cov2corr(iC),'parent',handles.ha(3))
axis(handles.ha(3),'square')
axis(handles.ha(3),'equal')
axis(handles.ha(3),'tight')
colorbar('peer',handles.ha(3))
colormap(cool)
VBA_title(handles.ha(3),'parameters'' correlation matrix')
set(handles.ha(3),'clim',[-1,1],'xdir','normal','ydir','reverse','xtick',[1:k],'ytick',[1:k])

% display contrast
handles.ha(4) = subplot(3,2,4,'parent',handles.hf,'visible','on');
switch type
    case 't'
        hp = bar(handles.ha(4),c,'facecolor',0.8*[1 1 1],'BarWidth',0.5);
        pos0 = get(handles.ha(4),'position');
        set(handles.ha(4),'position',[pos(1),pos0(2),pos(3),pos0(4)])
    case 'F'
        hi = imagesc(c','parent',handles.ha(4));
        colorbar('peer',handles.ha(4))
        set(handles.ha(4),'ytick',[1:1:m],'ylim',[0.5,m+0.5],'ygrid','on')
end
set(handles.ha(4),'xtick',[1:1:k],'xlim',[0.5,k+0.5],'ygrid','on')
xlabel(handles.ha(4),'independent variables')
VBA_title(handles.ha(4),'constrast')
box(handles.ha(4),'off')

% display test results (p-value and F-test)
pos1 = [0.4,0.03,0.6,0.02];
un = 'normalized';
handles.ht(1) = uicontrol('parent',handles.hf,'style','text','units',un,'position',pos1,'string',[],'fontsize',10,'backgroundcolor',[1 1 1],'HorizontalAlignment','left');

pos1 = pos1 + [0.2,0,0,0];
handles.ht(2) = uicontrol('parent',handles.hf,'style','text','units',un,'position',pos1,'string',strdf,'fontsize',10,'backgroundcolor',[1 1 1],'HorizontalAlignment','left');

% data selector
ud.type = type;
ud.all = all;
ud.y = y;
ud.X = X;
ud.yhat = yhat;
ud.pv = pv;
ud.stat = stat;
ud.Xnames = Xnames;
ud.Ynames = Ynames;
ud.type = type;
ud.handles = handles;
pos1 = [0.4,0.96,0.2,0.02];
handles.ht(2) = uicontrol('parent',handles.hf,'style','popupmenu','units',un,'position',pos1,'string',stryn,'fontsize',10,'backgroundcolor',0.8*[1 1 1],'HorizontalAlignment','left','userdata',ud,'callback',@myData);
feval(@myData,handles.ht(2),[])

% residuals inspection
str = 'inspect residuals';
pos2 = [0.15,0.032,0.2,0.02];
ud.handles = handles;
handles.ht(3) = uicontrol('parent',handles.hf,'style','pushbutton','units',un,'position',pos2,'string',str,'fontsize',10,'backgroundcolor',0.8*[1 1 1],'HorizontalAlignment','left','userdata',ud,'callback',@myResiduals);

all.handles = handles;

try,VBA_getSubplots ();end

function pBP = myBreuschPaganTest(e,X)
X = [X,ones(size(X,1),1)];
n = size(X,2);
c = [eye(n-1);zeros(1,n-1)];
pBP = GLM_contrast(X,e.^2,c,'F',0);

function myResiduals(e1,e2)
ud = get(e1,'userdata');
ind = get(ud.handles.ht(2),'value');
e = ud.y(:,ind) - ud.yhat(:,ind);
% derive residuals' empirical histogram
[n,x] = hist(e,16);
me = mean(e);
ve = var(e);
mie = min(e);
mae = max(e);
de = (mae-mie)*1e-3;
ge = mie:de:mae;
pe = exp(-0.5*(ge-me).^2./ve)./(sqrt(ve*2*pi));
dx = diff(x);
dx = dx(1);
[epe,ege] = VBA_empiricalDensity(e,1);
for i=1:length(x)
    in = find(ge>x(i)-dx/2&ge<x(i)+dx/2);
    sp(i) = mean(pe(in));
    in = find(ege>x(i)-dx/2&ege<x(i)+dx/2);
    sep(i) = mean(epe(in));
end
sc = sp*n'./(sp*sp');
esc = sep*n'./(sep*sep');
% look for heteroscedasticity
pBP = myBreuschPaganTest(e,ud.X);
for i=1:10
    ming = quantile(ud.yhat(:,ind),(i-1)*0.1);
    maxg = quantile(ud.yhat(:,ind),i*0.1);
    in = find(ud.yhat(:,ind)>=ming&ud.yhat(:,ind)<=maxg);
    myh(i) = mean(ud.yhat(in,ind));
    meh(i) = mean(e(in));
    seh(i) = std(e(in));
end
% display results
sn = ['GLM residuals: data #',num2str(ind)];
try
    sn = [sn,' (',ud.Ynames{ind},')'];
end
pos0 = get(ud.handles.hf,'position');
hf = figure('color',[1 1 1],'name',sn,'menubar','none');
pos = get(hf,'position');
set(hf,'position',[pos(1),pos(2),pos0(3),pos(4)])
ha = subplot(1,2,1,'parent',hf,'nextplot','add');
hb = bar(x,n,'parent',ha);
set(hb,'facecolor',0.8*[1 1 1])
plot(ha,ege,epe*esc,'k')
plot(ha,ge,pe*sc,'r')
xlabel(ha,'residual bins')
ylabel(ha,'# samples per bin')
legend(ha,{'histogram of residuals','empirical distribution of residuals','gaussian approximation'})
VBA_title(ha,'residuals'' normality')
xl = get(ha,'xlim');
yl = get(ha,'ylim');
ht = text(xl(1)+dx/2,yl(2),['Kolmogorov-Smirnov test: p=',num2str(ud.all.ks(ind),3)],'color','r','parent',ha);
ha = subplot(1,2,2,'parent',hf,'nextplot','add');
he = errorbar(myh,meh,seh,'parent',ha);
set(he,'color','k','marker','.')
mih = min(ud.yhat(:,ind));
mah = max(ud.yhat(:,ind));
plot(ha,[mih,mah],[me,me],'r')
plot(ha,[mih,mah],[me+sqrt(ve),me+sqrt(ve)],'r--')
plot(ha,[mih,mah],[me-sqrt(ve),me-sqrt(ve)],'r--')
yl = get(ha,'ylim');
ht = text(mih+(mah-mih)/32,yl(2),['Breusch-Pagan test: p=',num2str(pBP,3)],'color','r','parent',ha);
ht = text(mih+(mah-mih)/32,yl(2)-0.05*diff(yl),['E[residuals]=',num2str(me,3)],'color','r','parent',ha);
ht = text(mih+(mah-mih)/32,yl(2)-0.1*diff(yl),['V[residuals]=',num2str(ve,3)],'color','r','parent',ha);
set(ha,'xlim',[mih,mah])
xlabel(ha,'predicted data')
ylabel(ha,'residuals'' moments')
VBA_title(ha,'residuals'' homoscedasticity')
try;VBA_getSubplots ();end

function myData(e1,e2)
ud = get(e1,'userdata');
ind = get(e1,'value');

% plot data
cla(ud.handles.ha(1))
hp = plot(ud.handles.ha(1),ud.yhat(:,ind),ud.y(:,ind),'k.');
mi = min([ud.y(:,ind);ud.yhat(:,ind)]);
ma = max([ud.y(:,ind);ud.yhat(:,ind)]);
plot(ud.handles.ha(1),[mi,ma],[mi,ma],'r')
xx = 0.1*(ma-mi)+mi;
str = ['R^2=',sprintf('%2.3f',ud.all.R2(ind))];
if isequal(ud.type,'F')
    str = [str, ' [R^2[contrast]=',sprintf('%2.3f',ud.all.R2_a(ind)),']'];
end
text(xx,xx,str,'parent',ud.handles.ha(1),'color','r')
axis(ud.handles.ha(1),'tight')
xlabel(ud.handles.ha(1),'predicted data')
ylabel(ud.handles.ha(1),'observed data')
grid(ud.handles.ha(1),'on')
VBA_title(ud.handles.ha(1),'data alignement')

cla(ud.handles.ha(6))
hp = plot(ud.handles.ha(6),ud.yhat(:,ind),'linestyle','-','marker','.','color',[1,0,0]);
hp = plot(ud.handles.ha(6),ud.y(:,ind),'k.');
legend(ud.handles.ha(6),{'predicted data','observed data'})
axis(ud.handles.ha(6),'tight')
xlabel(ud.handles.ha(6),'data dimensions')
ylabel(ud.handles.ha(6),'data')
grid(ud.handles.ha(6),'on')
VBA_title(ud.handles.ha(6),'data fit')

% parameter estimates
cla(ud.handles.ha(2))
Vb = diag(ud.all.vhat(ind)*ud.all.iC);
k = size(ud.all.b,1);
for j=1:k
    hp = bar(ud.handles.ha(2),j,ud.all.b(j,ind),'facecolor',0.8*[1 1 1],'BarWidth',0.5);
    hcmenu = uicontextmenu;
    str = ['variable #',num2str(j)];
    if ~isempty(ud.Xnames)
        str = [str,' (',ud.Xnames{j},')'];
    end
    uimenu(hcmenu, 'Label',str);
    try % if removed for  computational burden...
        uimenu(hcmenu, 'Label',['p=',num2str(ud.all.pv(j,ind),'%3.3f')]);
        uimenu(hcmenu, 'Label',['F=',num2str(ud.all.stat(j,ind),'%3.3f')]);
        uimenu(hcmenu, 'Label',['dof=[',num2str(ud.all.df(j,1,ind)),',',num2str(ud.all.df(j,2,ind)),']']);
        R2 = VBA_FtoR2(ud.all.stat(j,ind),ud.all.df(j,1,ind),ud.all.df(j,2,ind));
        uimenu(hcmenu, 'Label',['R2=',num2str(R2*100,'%3.1f'),'%']);
    end
    set(get(hp,'children'),'uicontextmenu',hcmenu);
    set(hp,'uicontextmenu',hcmenu);
    hp = errorbar(ud.handles.ha(2),j,ud.all.b(j,ind),1.96*sqrt(Vb(j)),'r.'); %p=.05 confidence intervals
    %     hp = errorbar(ud.handles.ha(2),j,ud.all.b(j,ind),sqrt(Vb(j)),'r.');
    set(hp,'uicontextmenu',hcmenu);
end
set(ud.handles.ha(2),'xtick',[1:1:k],'xlim',[0.5,k+0.5],'ygrid','on')
xlabel(ud.handles.ha(2),'independent variables')
VBA_title(ud.handles.ha(2),'parameter estimates')
strp = ['p=',num2str(ud.pv(ind),'%3.3f'),' (',ud.type,'=',num2str(ud.stat(ind),'%3.3f'),')'];
set(ud.handles.ht(1),'string',strp);


