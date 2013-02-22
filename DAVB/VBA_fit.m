function fit = VBA_fit(posterior,out)
% derives standard model fit accuracy metrics
% function fit = VBA_fit(posterior,out)
% IN:
%   - posterior/out: output structures of VBA_NLStateSpaceModel.m
% OUT:
%   - fit: structure, containing the following fields:
%       .LL: log-likelihood of the model
%       .R2: coefficient of determination (the fraction of variance
%       unexplained is 1-R2)
%       .AIC: Akaike Information Criterion
%       .BIC: Bayesian Informaion Criterion

suffStat = out.suffStat;

% Log-likelihood
if ~out.options.binomial
    v = posterior.b_sigma./posterior.a_sigma;
    fit.LL = -0.5*out.suffStat.dy2/v;
    fit.ny = 0;
    for t=1:out.dim.n_t
        ldq = VBA_logDet(out.options.priors.iQy{t}/v);
        fit.ny = fit.ny + length(find(diag(out.options.priors.iQy{t})~=0));
        fit.LL = fit.LL - 0.5*ldq;
    end
    fit.LL = fit.LL - 0.5*fit.ny*log(2*pi);
else
    fit.LL = out.suffStat.logL;
    fit.ny = sum(1-out.options.isYout(:));
end

% coefficient of determination
SS_tot = sum((out.y(:)-mean(out.y(:))).^2);
SS_err = sum((out.y(:)-suffStat.gx(:)).^2);
fit.R2 = 1-(SS_err/SS_tot);

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
fit.BIC = fit.LL - 0.5*fit.ntot*log(fit.ny);


