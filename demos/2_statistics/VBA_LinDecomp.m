function [A,B,out,posterior] = VBA_LinDecomp(Y,k,options)
% linear decomposition of a matrix Y
% function [A,B,out] = VBA_LinDecomp(Y,k)
% This function operates a probabilistic linear decomposition of matrix Y,
% as follows: Y = A*B + cst + e, where A and B are unknow matrices, cst is
% a constant term and e is assumed to be random noise.
% IN:
%   - Y: nxp data matrix
%   - k: number of components (default=1)
%   - options: VBA options' structure (NB: priors are overwritten!)
% OUT:
%   - A: nxk loadings matrix
%   - B: kxp dual loadings matrix
%   - out/posterior: VBA output structures

try,k;catch,k=1;end
try,options;catch,options=[];end
[n,p] = size(Y);

% set parameter indices for model inversion
inG.n = k;
inG.size.X = [n,k];
inG.size.Y = [k,p];
for i=1:k
    inG.ind(i).X = (i-1)*n+1:i*n;
    inG.ind(i).Y = n*k+(i-1)*p+1:n*k+i*p;
end
inG.ind0 = max(inG.ind(k).Y)+1;
options.inG = inG;
g_fname = @g_LinDecomp;
dim.n = 0;
dim.n_theta = 0;
dim.n_phi = inG.ind0;
% set priors and invert model
options.priors.muPhi = 1e-3*ones(dim.n_phi,1);
options.priors.SigmaPhi = 1e2*eye(dim.n_phi);
[posterior,out] = VBA_hyperparameters(VBA_vec(Y),[],[],g_fname,dim,options);
% recover A and B matrices from VBA posterior pdf
A = NaN(n,k);
B = NaN(k,p);
for i=1:inG.n
    A(:,i) = posterior.muPhi(inG.ind(i).X);
    B(i,:) = posterior.muPhi(inG.ind(i).Y)';
end
% display results?
if out.options.DisplayWin
    Yh = A*B + posterior.muPhi(out.options.inG.ind0);
    name1 = 'observed data';
    name2 = 'fitted data';
    mi = min([Y(:);Yh(:)]);
    ma = max([Y(:);Yh(:)]);
    figname = ['VBA_LinDecomp (k=',num2str(k),'): R2=',num2str(out.fit.R2*100,3),'%'];
    hf = figure('color',[1 1 1],'menubar','none','name',figname);
    ha = subplot(3,2,1,'parent',hf);
    plot(ha,Y(:),Yh(:),'b.')
    grid(ha,'off')
    axis(ha,'tight')
    xlabel(ha,name1)
    ylabel(ha,name2)
    title(ha,[name1,' vs ',name2])
    box(ha,'off')
    ha = subplot(3,2,2,'parent',hf);
    imagesc(Y-Yh,'parent',ha),colorbar
    title(ha,[name1,' - ',name2])
    ha = subplot(3,2,3,'parent',hf);
    imagesc(Y,'parent',ha),colorbar
    title(ha,name1)
    set(ha,'clim',[mi,ma])
    ha = subplot(3,2,4,'parent',hf);
    imagesc(Yh,'parent',ha),colorbar
    title(ha,name2)
    set(ha,'clim',[mi,ma])
    SS_tot = sum((Y(:)-mean(Y(:))).^2);
    dR2 = NaN(k,1);
    yhi = NaN(n*p,k);
    for i=1:k
        ii = setdiff(1:k,i);
        Yhi = A(:,ii)*B(ii,:) + posterior.muPhi(out.options.inG.ind0);
        SS_err = sum((Y(:)-Yhi(:)).^2);
        R2i = 1-(SS_err/SS_tot);
        dR2(i) = out.fit.R2 - R2i;
        yhi(:,i) = VBA_vec(A(:,i)*B(i,:));
    end
    ha = subplot(3,2,5,'parent',hf);
    hb = bar(1:k,dR2);
    set(hb,'facecolor',0.8*[1 1 1])
    xlabel(ha,'component removed')
    ylabel(ha,'loss of explained variance')
    box(ha,'off')
    title(ha,'components'' contribution')
    set(ha,'xtick',1:k,'xlim',[0,k+1])
    if k >1
        ha = subplot(3,2,6,'parent',hf);
        C = corr(yhi);
        imagesc(C,'parent',ha)
        set(ha,'clim',[-1,1],'xtick',1:k,'ytick',1:k)
        colorbar('peer',ha)
        title(ha,'components'' correlation')
    end
    VBA_getSubplots ();
end



