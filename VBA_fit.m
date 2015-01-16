function fit = VBA_fit(posterior,out)
% derives standard model fit accuracy metrics
% function fit = VBA_fit(posterior,out)
% IN:
%   - posterior/out: output structures of VBA_NLStateSpaceModel.m
% OUT:
%   - fit: structure, containing the following fields:
%       .LL: log-likelihood of the model
%       .AIC: Akaike Information Criterion
%       .BIC: Bayesian Informaion Criterion
%       .R2: coefficient of determination (fraction of explained variance). 
%       .acc  balanced classification accuracy (fraction of correctly predicted outcomes).

suffStat = out.suffStat;


gsi = find([out.options.sources.type]==0);
for i=1:length(gsi)
    
    si=gsi(i);
    
    % Log-likelihood
    v(i) = posterior.b_sigma(i)/posterior.a_sigma(i);
    fit.LL(si) = -0.5*out.suffStat.dy2(i)/v(i);
    fit.ny(si) = 0;
    for t=1:out.dim.n_t
        ldq = VBA_logDet(out.options.priors.iQy{t,i}/v(i));
        fit.ny(si) = fit.ny(si) + length(find(diag(out.options.priors.iQy{t,i})~=0));
        fit.LL(si) = fit.LL(si) + 0.5*ldq;
    end
    fit.LL(si) = fit.LL(si) - 0.5*fit.ny(si)*log(2*pi);    
    
    % coefficient of determination
%     if isfield(out.options,'sources')

        idx = out.options.sources(si).out;
        y_temp = out.y(idx,:);
        y_temp = y_temp(out.options.isYout(idx,:) == 0);
        
        gx_temp = suffStat.gx(idx,:);
        gx_temp = gx_temp(out.options.isYout(idx,:) == 0);
        
        SS_tot = sum((vec(y_temp)-mean(vec(y_temp))).^2);
        SS_err = sum((vec(y_temp)-vec(gx_temp)).^2);
        fit.R2(si) = 1-(SS_err/SS_tot);    
%     end
        fit.acc(si) = NaN;
end



bsi = find([out.options.sources.type]==1);
for i=1:length(bsi)
    si=bsi(i);
    idx = out.options.sources(si).out;

    fit.LL(si) = out.suffStat.logL(si);
    fit.ny(si) = sum(1-vec(out.options.isYout(idx,:)));
    
    % balanced accuracy
        y_temp = out.y(idx,:);
        y_temp = y_temp(out.options.isYout(idx,:) == 0);
        
        gx_temp = suffStat.gx(idx,:);
        gx_temp = gx_temp(out.options.isYout(idx,:) == 0);
        
        bg = gx_temp>.5; % binarized model predictions
        tp = sum(vec(y_temp).*vec(bg)); % true positives
        fp = sum(vec(1-y_temp).*vec(bg)); % false positives
        fn = sum(vec(y_temp).*vec(1-bg)); % false positives
        tn = sum(vec(1-y_temp).*vec(1-bg)); %true negatives
        P = tp + fn;
        N = tn + fp;
        fit.R2(si) = (tp+tn)./(P+N);
        fit.acc(si) = balanced_accuracy(suffStat.gx(idx,:),out.y(idx,:),out.options.isYout(idx,:));
    
end

msi = find([out.options.sources.type]==2);
for i=1:length(msi)
    
    si=msi(i);
    fit.LL(si) = out.suffStat.logL(si);
    
    idx = out.options.sources(si).out;

    fit.ny(si) = sum(1-any(out.options.isYout(idx,:)));
   
    
    fit.R2(si) = NaN;
    fit.acc(si) = multinomial_accuracy(suffStat.gx(idx,:),out.y(idx,:),out.options.isYout(idx,:));
    
end

% AIC/BIC
fit.ntot = 0;
if out.dim.n_phi > 0
    indIn = out.options.params2update.phi;
    fit.ntot = fit.ntot + length(indIn);
end
if out.dim.n_theta > 0
    indIn = out.options.params2update.theta;
    fit.ntot = fit.ntot + length(indIn);
end
if out.dim.n > 0  && ~isinf(out.options.priors.a_alpha) && ~isequal(out.options.priors.b_alpha,0)
    for t=1:out.dim.n_t
        indIn = out.options.params2update.x{t};
        fit.ntot = fit.ntot + length(indIn);
    end
end
fit.AIC = sum(fit.LL) - fit.ntot;
fit.BIC = sum(fit.LL) - 0.5*fit.ntot.*log(sum(fit.ny));

