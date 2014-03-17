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
%       .R2: if data is continuous, R2 = coefficient of determination
%       (fraction of explained variance). If data is binary, R2 = balanced
%       classification accuracy (fraction of correctly predicted outcomes).

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
        SS_tot = sum((vec(out.y(idx,:))-mean(vec(out.y(idx,:)))).^2);
        SS_err = sum((vec(out.y(idx,:))-vec(suffStat.gx(idx,:))).^2);
        fit.R2(si) = 1-(SS_err/SS_tot);
%     end
end



bsi = find([out.options.sources.type]~=0);
for i=1:length(bsi)
    si=bsi(i);
    fit.LL(si) = out.suffStat.logL(si);
    fit.ny(si) = sum(1-out.options.isYout(:));
    
    % balanced accuracy
%     if isfield(out.options,'sources')
        idx = out.options.sources(si).out;
        bg = out.suffStat.gx(idx,:)>.5; % binarized model predictions
        tp = sum(vec(out.y(idx,:)).*vec(bg)); % true positives
        fp = sum(vec(1-out.y(idx,:)).*vec(bg)); % false positives
        fn = sum(vec(out.y(idx,:)).*vec(1-bg)); % false positives
        tn = sum(vec(1-out.y(idx,:)).*vec(1-bg)); %true negatives
        P = tp + fn;
        N = tn + fp;
        fit.R2(si) = 0.5*(tp./P + tn./N);
        fit.acc(si) = (tp+tn)./(P+N);
%     end
    
end

% if ~isfield(out.options,'sources')
%     if ~out.options.binomial
%         % coefficient of determination
%         SS_tot = sum((out.y(:)-mean(out.y(:))).^2);
%         SS_err = sum((out.y(:)-suffStat.gx(:)).^2);
%         fit.R2 = 1-(SS_err/SS_tot);
%     else
%         % balanced accuracy
%         bg = out.suffStat.gx>.5; % binarized model predictions
%         tp = sum(vec(out.y).*vec(bg)); % true positives
%         fp = sum(vec(1-out.y).*vec(bg)); % false positives
%         fn = sum(vec(out.y).*vec(1-bg)); % false positives
%         tn = sum(vec(1-out.y).*vec(1-bg)); %true negatives
%         P = tp + fn;
%         N = tn + fp;
%         fit.R2 = 0.5*(tp./P + tn./N);
%         fit.acc = (tp+tn)./(P+N);
%     end
% end
%     


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
fit.AIC = fit.LL - fit.ntot;
fit.BIC = fit.LL - 0.5*fit.ntot.*log(fit.ny);


