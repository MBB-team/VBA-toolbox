% demo for binary data classification

close all
clear variables

dim.p = 32;
dim.n_t = 1;
dim.n_phi = 16;
dim.n_theta = 0;
dim.n = 0;


g_fname = @g_classif;
options.sources = struct('type',1 ,'out', 1); % one binomial observation;
options.priors.muPhi = zeros(dim.n_phi,1);
options.priors.SigmaPhi = 1e0*eye(dim.n_phi);
options.isYout = zeros(dim.p,1);
options.DisplayWin = 0;
options.verbose = 0;

Nmcmc = 64;
p = cell(2,2,Nmcmc);
o = cell(2,2,Nmcmc);
F = zeros(2,2,Nmcmc);
ner = zeros(2,2,Nmcmc); % proportion of correct predictions
mner = zeros(2,Nmcmc); % maximum performance rate

fprintf(1,['MCMC simulations... ']);
fprintf(1,'%6.2f %%',0)

for ii=1:Nmcmc
    options.inG.X = randn(dim.n_phi-1,dim.p);
    for i=1:2
        % simulate data with and without real mapping
        phi = (2-i)*randn(dim.n_phi,1);
        [y,x,x0,eta,e] = VBA_simulate (dim.n_t,[],g_fname,[],phi,[],[],[],options);
        g = y-e;
        g = g>0.5; % denoised data
        mner(i,ii) = sum(g.*y + (1-g).*(1-y))./dim.p; % max performance rate
        for j=1:2
            % invert model with and without the 2nd half of the data
            options.isYout(dim.p/2:dim.p) = 2-j;
            [p{j,i,ii},o{j,i,ii}] = VBA_NLStateSpaceModel(y,[],[],g_fname,dim,options);
            F(j,i,ii) = o{j,i,ii}.F - VBA_LMEH0(y,o{j,i,ii}.options);
            % proportion of correct predictions on 2nd half of the data
            gx = o{j,i,ii}.suffStat.gx(dim.p/2:dim.p);
            g0 = y(dim.p/2:dim.p);
            ner(j,i,ii) = sum(gx.*g0 + (1-gx).*(1-g0))./(dim.p/2);
        end
        
    end
    
    fprintf(1,repmat('\b',1,8))
    fprintf(1,'%6.2f %%',100*ii/Nmcmc)
    
    
end

fprintf(1,repmat('\b',1,8))
fprintf(' OK.')
fprintf('\n')


hf = figure('color',[1 1 1]);
pos = get(hf,'position');
set(hf,'position',pos.*[1 1 1.5 1]);
test = {'test','train'};
data = {'half dataset','full dataset'};
model = {'H1','H0'};
ylim = [0;0];

for j=1:2
    
    ha(j,1) = subplot(2,2,(j-1)*2+1,'parent',hf);
    mr = squeeze(mean(ner(j,:,:),3));
    vr = squeeze(var(ner(j,:,:),[],2)./Nmcmc);
    plotUncertainTimeSeries(mr',vr',[],ha(j,1));
    set(ha(j,1),'xlim',[0,3],'xtick',[1,2],'xticklabels',model)
    xlabel(ha(j,1),'type of simulated data')
    ylabel(ha(j,1),['P[correct prediction]'])
    title(ha(j,1),test{j})
    box(ha(j,1),'off')
    hold(ha(j,1),'on')
    plot(ha(j,1),[0,3],[0.5,0.5],'r--')
    mmr = mean(mner,2);
    plot(ha(j,1),[1,2],mmr,'go')
    
    ha(j,2) = subplot(2,2,(j-1)*2+2,'parent',hf);
    dF = F(j,:,:);
    mdF = squeeze(mean(dF,3));
    vdF = squeeze(var(dF,[],3)./Nmcmc);
    plotUncertainTimeSeries(mdF',vdF',[],ha(j,2));
    set(ha(j,2),'xlim',[0,3],'xtick',[1,2],'xticklabels',model)
    xlabel(ha(j,2),'type of simulated data')
    ylabel(ha(j,2),['log p(y|H1) - log p(y|H0)'])
    title(ha(j,2),'evidence for a mapping')
    box(ha(j,2),'off')
    title(ha(j,2),data{j})
    
    ylim(1) = min([ylim(1),get(ha(j,2),'ylim')]);
    ylim(2) = max([ylim(2),get(ha(j,2),'ylim')]);
end

for j=1:2
    set(ha(j,2),'ylim',ylim)
end


VBA_getSubplots ();



