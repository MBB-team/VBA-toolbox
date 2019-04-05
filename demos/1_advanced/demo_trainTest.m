% demo for train/test cross-validation schemes

close all
clear all

dim.p = 32;
dim.n_t = 1;
dim.n_phi = 4;
dim.n_theta = 0;
dim.n = 0;


g_fname = @g_classif;
options.sources = struct('type',1 ,'out', 1); % one binomial observation;
options.priors.muPhi = zeros(dim.n_phi,1);
options.priors.SigmaPhi = 1e0*eye(dim.n_phi);
options.DisplayWin = 0;
options.verbose = 0;

Nmcmc = 64;

train = zeros(4,Nmcmc);
test = zeros(4,Nmcmc);
train_b = zeros(4,Nmcmc);
test_b = zeros(4,Nmcmc);

folds = cell(4,1); %4-folds train/test paritions
for i=1:4
    folds{i} = (i-1)*dim.p/4+1:i*dim.p/4;
end

fprintf(1,['MCMC simulations... ']);
fprintf(1,'%6.2f %%',0)


for ii=1:Nmcmc
    % pick a model at random (as well as its parameters)
    options.inG.X = randn(dim.n_phi-1,dim.p);
    phi = randn(dim.n_phi,1);
    [y,x,x0,eta,e] = VBA_simulate (dim.n_t,[],g_fname,[],phi,[],[],[],options,[]);
    g = y-e;
    g = g>0.5; % denoised data
    % 4-fold train/test
    for j=1:4
        options.isYout = zeros(dim.p,1);
        options.isYout(folds{j}) = 1;
        [p,o] = VBA_NLStateSpaceModel(y,[],[],g_fname,dim,options);
        
        % measure train prediction accuracies
        gx0 = o.suffStat.gx(setdiff(1:dim.p,folds{j}));
        y0 = y(setdiff(1:dim.p,folds{j}));
        train(j,ii) = sum(gx0.*y0 + (1-gx0).*(1-y0))./(length(y0));
        train_b(j,ii) = 0;
        if length(find(y0==1)) >0
            train_b(j,ii) = train_b(j,ii)+ 0.5*sum(gx0.*y0)./length(find(y0==1));
        end
        if length(find(y0==0)) >0
            train_b(j,ii) = train_b(j,ii)+ 0.5*sum((1-gx0).*(1-y0))./length(find(y0==0));
        end
        
        % measure test prediction accuracies
        gx0 = o.suffStat.gx(folds{j});
        y0 = y(folds{j});
        test(j,ii) = sum(gx0.*y0 + (1-gx0).*(1-y0))./(length(y0));
        test_b(j,ii) = 0;
        if length(find(y0==1)) >0
            test_b(j,ii) = test_b(j,ii)+ 0.5*sum(gx0.*y0)./length(find(y0==1));
        end
        if length(find(y0==0)) >0
            test_b(j,ii) = test_b(j,ii)+ 0.5*sum((1-gx0).*(1-y0))./length(find(y0==0));
        end
        
    end
    
    fprintf(1,repmat('\b',1,8))
    fprintf(1,'%6.2f %%',100*ii/Nmcmc)
    
end

fprintf(1,repmat('\b',1,8))
fprintf(' OK.')
fprintf('\n')


hf = figure('color',[1 1 1]);

ha = subplot(3,2,1,'parent',hf);
mt = mean(train,1);
plotUncertainTimeSeries(mean(mt),var(mt),1,ha)
hold(ha,'on')
mt = mean(test,1);
plotUncertainTimeSeries(mean(mt),var(mt),2,ha)
set(ha,'xtick',[1,2],'xticklabels',{'train','test'})
title(ha,'raw accuracy')

ha = subplot(3,2,2,'parent',hf);
mt = mean(train_b,1);
plotUncertainTimeSeries(mean(mt),var(mt),1,ha)
hold(ha,'on')
mt = mean(test_b,1);
plotUncertainTimeSeries(mean(mt),var(mt),2,ha)
set(ha,'xtick',[1,2],'xticklabels',{'train','test'})
title(ha,'balanced accuracy')

ha = subplot(3,2,3,'parent',hf);
md = mean(train,1) - mean(test,1);
[n,x] = hist(md);
bar(x,n,'parent',ha)
title(ha,'E[acc(train)]-E[acc(test)]')

ha = subplot(3,2,4,'parent',hf);
mdb = mean(train_b,1) - mean(test_b,1);
[n,x] = hist(mdb);
bar(x,n,'parent',ha)
title(ha,'E[acc(train)]-E[acc(test)] (balanced)')

ha = subplot(3,2,5,'parent',hf);
mv = var(train,[],1) - var(test,[],1);
[n,x] = hist(mv);
bar(x,n,'parent',ha)
title(ha,'V[acc(train)]-V[acc(test)]')

ha = subplot(3,2,6,'parent',hf);
mvb = var(train_b,[],1) - var(test_b,[],1);
[n,x] = hist(mvb);
bar(x,n,'parent',ha)
title(ha,'V[acc(train)]-V[acc(test)] (balanced)')


VBA_getSubplots ();



