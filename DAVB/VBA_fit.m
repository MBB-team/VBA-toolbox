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
gsi = find([out.options.sources.type]==0);
for i=1:length(gsi)
    si=gsi(i);
    v(i) = posterior.b_sigma(i)/posterior.a_sigma(i);
    fit.LL(si) = -0.5*out.suffStat.dy2(i)/v(i);
    fit.ny(si) = 0;
    for t=1:out.dim.n_t
        ldq = VBA_logDet(out.options.priors.iQy{t,i}/v(i));
        fit.ny(si) = fit.ny(si) + length(find(diag(out.options.priors.iQy{t,i})~=0));
        fit.LL(si) = fit.LL(si) - 0.5*ldq;
    end
    fit.LL(si) = fit.LL(si) - 0.5*fit.ny(si)*log(2*pi);
end
bsi = find([out.options.sources.type]~=0);
for i=1:length(bsi)
    si=bsi(i);
    fit.LL(si) = out.suffStat.logL(si);
    fit.ny(si) = sum(1-out.options.isYout(:));
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
fit.BIC = fit.LL - 0.5*fit.ntot.*log(fit.ny);


