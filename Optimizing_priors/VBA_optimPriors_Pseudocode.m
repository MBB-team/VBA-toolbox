function [priors] = VBA_optimPriors_Pseudocode(M,u,partition,density)


% E[e(y)] using either MCMC or Laplace

% if MCMC:


% loop over priors----------------------------------------

for k = 1:lngth(gridp)
    
    for m=1:length(partition)
        for j=1:length(partition{m})
            ind = partition{m}(j);
            opt = M{ind}.options;
            opt.priors = density;
            [pX,gX,pY,gY,X,Y] = get_MCMC_predictiveDensity(M{ind}.f_fname,M{ind}.g_fname,u,n_t,opt,M{ind}.dim,N);
            for i=1:N
                for m'=1:length(M)
                    opt = M{m'}.options;
                    opt.priors = Priors(priors,k,gridp)
                    [posterior,out] = VBA_NLStateSpaceModel(Y(:,:,i),u,M{m'}.f_fname,M{m'}.g_fname,M{m'}.dim,opt)
                    F(j,m',i) = out.F
                end
                e(i) = scoreError(F(j,:,i),partition,m);
            end
            Ee(m,j) = mean(e);
        end
    end
    em(k) = mean(Ee(:));
    
    
%     for m=1:length(partition)
%         for j=1:length(partition{m})
%             ind = partition{m}(j);
%             opt = M{ind}.options;
%             opt.priors = density;
%             [muy,Vy] = VBA_getLaplace(u,M{ind}.f_fname,M{ind}.g_fname,M{ind}.dim,opt)
%             em = getEfromY(muy,u,M,m,k,gridp);
%             dedy = numericDiff(@getFfromY,1,muy,u,M,m);
%             Ee(m,j) = em + 0.5*dedy'*Vy*dedy;
%         end
%     end
%     em(k) = mean(Ee);
%     
    
end

% loop over priors----------------------------------------

% find priors that minimizes em



function e = getEfromY(y,u,M,m,k,gridp)
for m'=1:length(M)
    opt = M{m'}.options;
    opt.priors = Priors(priors,k,gridp)
    [posterior,out] = VBA_NLStateSpaceModel(y,u,M{m'}.f_fname,M{m'}.g_fname,M{m'}.dim,opt)
    F(m') = out.F
end
e = scoreError(F,partition,m);

function e = scoreError(F,partition,true);
F = F - max(F);
p = exp(F)./sum(exp(F));
for i=1:length(partition)
    pp(i) = sum(p(partition{i}));
end
e = 1 - pp(true);


function priors = Priors(priors,k,gridp)
