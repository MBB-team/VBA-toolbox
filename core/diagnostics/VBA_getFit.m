function fit = VBA_getFit(posterior,out)
% derives standard model fit accuracy metrics
% function fit = VBA_getFit(posterior,out)
% IN:
%   - posterior/out: output structures of VBA_NLStateSpaceModel.m
% OUT:
%   - fit: structure, containing the following fields:
%       .LL: log-likelihood of the model
%       .AIC: Akaike Information Criterion
%       .BIC: Bayesian Informaion Criterion
%       .R2: coefficient of determination (fraction of explained variance). 
%       .acc: classification accuracy (fraction of correctly predicted outcomes).
%       .bacc: balanced classification accuracy.
%       .ny: effective sample size (total data dimension - #excluded data points)
%       .np: effective number of unknown model variables (total #params - #fixed params)
%
% [Note]: there was a change in this function on 09-06-2017. Before, the
% field .R2 used to report the coef of determination for continuous data
% and the balanced accuracy for binary data...

suffStat = out.suffStat;

% 0- effective number of unknown model variables
fit.np = 0;
if out.dim.n_phi > 0
    indIn = out.options.params2update.phi;
    fit.np = fit.np + length(indIn);
end
if out.dim.n_theta > 0
    indIn = out.options.params2update.theta;
    fit.np = fit.np + length(indIn);
end
if out.dim.n > 0
    indIn = out.options.params2update.x0;
    fit.np = fit.np + length(indIn);
    if ~isinf(out.options.priors.a_alpha) && ~isequal(out.options.priors.b_alpha,0)
        for t=1:out.dim.n_t
            indIn = out.options.params2update.x{t};
            fit.np = fit.np + length(indIn);
        end
    end
end

% 1- gaussian sources: goodness-of-fit
gsi = find([out.options.sources.type]==0);
for i=1:length(gsi)
    si=gsi(i);
    idx = out.options.sources(si).out;
    % sample size
    fit.ny(si) = sum(1-vec(out.options.isYout(idx,:)));
    % log-likelihood
    if out.options.UNL % to be rationalized...
        fit.LL = out.suffStat.logL;
    else
        v(i) = posterior.b_sigma(i)/posterior.a_sigma(i);
        fit.LL(si) = -0.5*out.suffStat.dy2(i)/v(i);
        for t=1:out.dim.n_t
            ldq = VBA_logDet(out.options.priors.iQy{t,i}/v(i));
            fit.LL(si) = fit.LL(si) + 0.5*ldq;
        end
        fit.LL(si) = fit.LL(si) - 0.5*fit.ny(si)*log(2*pi);
    end
    % AIC/BIC
    fit.AIC(si) = fit.LL(si) - fit.np;
    fit.BIC(si) = fit.LL(si) - 0.5*fit.np.*log(fit.ny(si));
    % coefficient of determination
    y_temp = out.y(idx,:);
    y_temp = y_temp(out.options.isYout(idx,:) == 0);
    gx_temp = suffStat.gx(idx,:);
    gx_temp = gx_temp(out.options.isYout(idx,:) == 0);
    SS_tot = sum((vec(y_temp)-mean(vec(y_temp))).^2);
    SS_err = sum((vec(y_temp)-vec(gx_temp)).^2);
    fit.R2(si) = 1-(SS_err/SS_tot);
    % classification accuracies [irrelevant]
    fit.acc(si) = NaN;
    fit.bacc(si) = NaN;
end


% 2- binomial sources: goodness-of-fit
bsi = find([out.options.sources.type]==1);
for i=1:length(bsi)
    si=bsi(i);
    idx = out.options.sources(si).out;
    % sample size
    fit.ny(si) = sum(1-vec(out.options.isYout(idx,:)));
    % log-likelihood
    fit.LL(si) = out.suffStat.logL(si);
    % AIC/BIC
    fit.AIC(si) = fit.LL(si) - fit.np;
    fit.BIC(si) = fit.LL(si) - 0.5*fit.np.*log(fit.ny(si));
    % coefficient of determination
    y_temp = out.y(idx,:);
    y_temp = y_temp(out.options.isYout(idx,:) == 0);
    gx_temp = suffStat.gx(idx,:);
    gx_temp = gx_temp(out.options.isYout(idx,:) == 0);
    SS_tot = sum((vec(y_temp)-mean(vec(y_temp))).^2);
    SS_err = sum((vec(y_temp)-vec(gx_temp)).^2);
    fit.R2(si) = 1-(SS_err/SS_tot);
    % classification accuracies    
    [fit.acc(si), fit.bacc(si)] = VBA_accuracy (suffStat.gx(idx,:), out.y(idx,:), 1, out.options.isYout(idx,:));
    
end

% 3- multinomial sources: goodness-of-fit
msi = find([out.options.sources.type]==2);
for i=1:length(msi)
    si=msi(i);
    idx = out.options.sources(si).out;
    % sample size
    fit.ny(si) = sum(1-any(out.options.isYout(idx,:)));
    % log-likelihood
    fit.LL(si) = out.suffStat.logL(si);
    % AIC/BIC
    fit.AIC(si) = fit.LL(si) - fit.np;
    fit.BIC(si) = fit.LL(si) - 0.5*fit.np.*log(fit.ny(si));
    % coefficient of determination
    y_temp = out.y(idx,:);
    y_temp = y_temp(out.options.isYout(idx,:) == 0);
    gx_temp = suffStat.gx(idx,:);
    gx_temp = gx_temp(out.options.isYout(idx,:) == 0);
    SS_tot = sum((vec(y_temp)-mean(vec(y_temp))).^2);
    SS_err = sum((vec(y_temp)-vec(gx_temp)).^2);
    fit.R2(si) = 1-(SS_err/SS_tot);
    % classification accuracies [to be rationalized!]
    [fit.acc(si), fit.bacc(si)] = VBA_accuracy (suffStat.gx(idx,:), out.y(idx,:), 2, out.options.isYout(idx,:));
    
end






