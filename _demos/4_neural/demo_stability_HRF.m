% explore balloon model stability

clear all
close all
clc

xg1 = [-5:0.2:5];
xg2 = xg1;
n = length(xg1);
lm = zeros(n,n);
kas = [-4:0.2:4];
kaf = [-3:0.2:3];

try

    load demo_stability_hrf.mat
    % get max eigenvalue
    P = [0;0;0;0;0;0];
    for i=1:n
        for j=1:n
            J = VBA_numericDiff(@f_HRF,1,[xg1(i);xg2(j);0;0],P,0,[]);
            lm(i,j) = max(real(eig(J)));
        end
    end

    % get separatrix
    stab = zeros(n,n);
    stab(lm>=0) = -1;
    L = del2(stab);
    sep = zeros(n,n);
    sep(L<0)=1;
    p = length(find(L==1));
    y = [];
    for i=1:n
        for j=1:n
            if sep(i,j)==1
                y = [y,[xg1(i);xg2(j)]];
            end
        end
    end
    g_fname = @g_exp;
    options.inG.x = y(2,:)';
    options.priors.muPhi = [1;1;0];
    options.priors.SigmaPhi = eye(3);
    options.priors.SigmaPhi(1,1) = 0;
    options.priors.SigmaPhi(3,3) = 0;
    options.verbose = 0;
    options.DisplayWin = 0;
    dim.n_theta         = 0;
    dim.n_phi           = 3;
    dim.n               = 0;
    [posterior,out] = VBA_NLStateSpaceModel(y(1,:)',[],[],g_fname,dim,options);

catch

    % loop over different kas/kaf
    phi = zeros(length(kaf),length(kas));
    for k=1:length(kas)
        for l=1:length(kaf)
            k,l
            % get max eigenvalue
            P = [0;0;kaf(l);kas(k);0;0];
            for i=1:n
                for j=1:n
                    J = VBA_numericDiff(@f_HRF,1,[xg1(i);xg2(j);0;0],P,0,[]);
                    lm(i,j) = max(real(eig(J)));
                end
            end

            % get separatrix
            stab = zeros(n,n);
            stab(lm>=0) = -1;
            L = del2(stab);
            sep = zeros(n,n);
            sep(L<0)=1;
            p = length(find(L==1));
            y = [];
            for i=1:n
                for j=1:n
                    if sep(i,j)==1
                        y = [y,[xg1(i);xg2(j)]];
                    end
                end
            end
            g_fname = @g_exp;
            options.inG.x = y(2,:)';
            options.priors.muPhi = [1;1;0];
            options.priors.SigmaPhi = eye(3);
            options.priors.SigmaPhi(1,1) = 0;
            options.priors.SigmaPhi(3,3) = 0;
            options.verbose = 0;
            options.DisplayWin = 0;
            dim.n_theta         = 0;
            dim.n_phi           = 3;
            dim.n               = 0;
            [posterior,out] = VBA_NLStateSpaceModel(y(1,:)',[],[],g_fname,dim,options);
            phi(l,k) = posterior.muPhi(2);

        end
    end

end

% store for display purposes
x = options.inG.x;
gx = out.suffStat.gx;
PHI = posterior.muPhi;

gkaf = kron(ones(1,length(kas)),kaf)';
gkaf = gkaf - mean(gkaf);
gkas = kron(kas,ones(1,length(kaf)))';
gkas = gkas - mean(gkas);
options = [];
g_fname = @g_exp2d;
options.inG.gkaf = gkaf;
options.inG.gkas = gkas;
options.priors.muPhi = ones(1,1);
options.priors.SigmaPhi = 1e4*eye(1);
options.verbose = 1;
options.DisplayWin = 1;
dim.n_theta         = 0;
dim.n_phi           = 1;
dim.n               = 0;
[posterior,out] = VBA_NLStateSpaceModel(...
    phi(:),[],[],g_fname,dim,options);
I = reshape(out.suffStat.gx,length(kaf),length(kas));


% display results

hf = figure('color',[1 1 1],'position',[649,40,767,1105],'menubar','none');
ha = subplot(3,2,1,'parent',hf);
hi = imagesc(flipud(lm),'parent',ha);
xlabel(ha,'x2: blood inflow')
ylabel(ha,'x1: vasodilatory signal')
title(ha,'max eigenvalue')
colorbar('peer',ha);

ha = subplot(3,2,2,'parent',hf);
hi = plot(y(2,:),y(1,:),'.','parent',ha);
xlabel(ha,'x2: blood inflow')
ylabel(ha,'x1: vasodilatory signal')
set(ha,...
    'xlim',[min(xg1),max(xg1)],...
    'ylim',[min(xg1),max(xg1)],...
    'nextplot','add')
axis(ha,'square')
grid(ha,'on')
title(ha,...
    ['separatrix: kas=',...
    num2str(P(4)),...
    ' ; kaf=',...
    num2str(P(3))])
set(ha,'nextplot','add')
[xs,is] = sort(x);
ys = gx(is);
plot(ha,xs,ys,'k','parent',ha)
plot(0,0,'r+','parent',ha);
str = ['x(1) = ',...
    'theta',...
    '*exp(',...
    'x(2)) : theta = ',...
    num2str(PHI(2))];
legend(ha,{'numerical',str,'equilibrium'})

ha = subplot(3,2,3,'parent',hf);
hi = imagesc(phi,'parent',ha);
colorbar('peer',ha);
xlabel(ha,'kas')
ylabel(ha,'kaf')
title('separatrix parameter (theta)')

ha = subplot(3,2,4,'parent',hf);
hi = plot(kas,phi','parent',ha);
xlabel(ha,'kas')
ylabel(ha,'separatrix parameter (theta)')
grid(ha,'on');
axis(ha,'tight')
for i=1:length(kaf)
    leg{i} = ['kaf=',num2str(kaf(i))];
end
legend(ha,leg)
title(ha,'separatrix parameter (theta)')

ha = subplot(3,2,5,'parent',hf);
hi = imagesc(I,'parent',ha);
colorbar('peer',ha);
xlabel(ha,'kas')
ylabel(ha,'kaf')
title('fitted separatrix parameter (theta)')

ha = subplot(3,2,6,'parent',hf);
hi = plot(phi(:),I(:),'k.','parent',ha);
mi = min([phi(:);I(:)]);
ma = max([phi(:);I(:)]);
set(ha,'nextplot','add')
plot(ha,[mi ma],[mi ma],'r')
xlabel(ha,'sep param')
ylabel(ha,'fitted sep param (theta)')
title('sep param model fit')

str = ['theta = ',...
    num2str(posterior.muPhi(1)),...
    '*exp(kaf/2 - |kas-kaf/2|)'];
legend(ha,{'numerical',str})
grid(ha,'on');
set(ha,'xlim',[mi ma],'ylim',[mi ma])
