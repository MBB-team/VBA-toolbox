function [gx] = fitFixedForm(y,u,ue,verbose)
% fit and extrapolate monotonic time series using fixed form models
% function [gx] = fitFixedForm(y,u,ue)
% IN:
%   - y: the raw time series
%   - u: the time grid
%   - ue: the time grid, on which the time series is to be extrapolated
% OUT:
%   - gx: extrapolated time series


try
    verbose;
catch
    verbose = 1;
end

g_fnames = {'cst','lin','qua','pow','cvg','cvg2'};
inG.diff = 0;
inG.div = 0;
inG.int = 0;
inG.log = 0;
inG.u = u;
options.inG = inG;
options.verbose = 0;
options.DisplayWin = verbose;
priors.a_sigma = 1e2;
priors.b_sigma = 1e0;
priors.muPhi = zeros(3,1);
priors.SigmaPhi = 1e2*eye(3);
priors.SigmaPhi(1,1) = 0;
options.priors      = priors;
dim.n_theta         = 0;
dim.n_phi           = 3;
dim.n               = 0;

nmods = length(g_fnames);

F = zeros(nmods,1);
F(3) = -Inf;
n = size(y,1);
yh = zeros(n,nmods);
ay = augmentgx(y,options.inG);
ay = ay(:);
for i=setdiff(1:nmods,3)
    disp(['model ',g_fnames{i},' (',num2str(i),'/',num2str(nmods),')...'])
    [p{i},o{i}] = VBA_NLStateSpaceModel(ay,u,[],g_fnames{i},dim,options);
    F(i) = o{i}.F;
    yh(:,i) = o{i}.suffStat.gx(1:n);
end
[mF,iF] = max(F);

% gx = yh(:,iF);
P = p{iF}.muPhi;
gx = feval(g_fnames{iF},[],P,ue,options.inG);

if ~verbose
    return
end


hf = figure('color',[1 1 1]);
ha = subplot(2,2,1,'parent',hf,'nextplot','add');
hb = bar(F,'parent',ha);
plot(iF,mF,'ro')
set(ha,'xtick',1:nmods,'xticklabels',g_fnames)

ha(2) = subplot(2,2,2,'parent',hf,'nextplot','add');
plot(ha(2),u,y,'k--')
plot(ha(2),u,yh)
legend(ha(2),cat(2,'data',g_fnames))


if iF == 1
    str = ['g(x)=', num2str(P(1),3)];
elseif iF == 2
    str = ['g(x)=', num2str(P(1),3),' + ',num2str(P(2),3),'x'];
elseif iF == 3
    str = ['g(x)=', num2str(P(1),3),' + ',num2str(P(2),3),'x',' + ',num2str(P(3),3),'x^2'];
elseif iF == 4
    % gx = P(1) + P(2).*u.^P(3);
    str = ['g(x)=', num2str(P(1),3),' + ',num2str(P(2),3),'x.^',num2str(P(3),3)];
elseif iF == 5
    % gx = P(1) + P(2)*(1-exp(-exp(P(3))*u));
    str = ['g(x)=', num2str(P(1),3),' + ',num2str(P(2),3),'(1-exp(-',num2str(exp(P(3)),3),' x))'];
elseif iF == 6
    % gx = P(1) + P(2)*exp(-exp(P(3))*u);
    str = ['g(x)=', num2str(P(1),3),' + ',num2str(P(2),3),'exp(-',num2str(exp(P(3)),3),' x)'];
end

set(hf,'name',[g_fnames{iF},': ',str]);

ay = augmentgx(y,inG);
agx = augmentgx(yh(:,iF),inG);
ha(3) = subplot(2,2,3,'parent',hf,'nextplot','add');
plot(ha(3),u,ay,'k--')
plot(ha(3),u,agx)
getSubplots




