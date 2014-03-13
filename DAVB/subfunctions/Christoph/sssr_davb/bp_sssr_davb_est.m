% Fits the Belief-Precision model to observed data. Top-level function.
% 
% Usage:
%     est = bp_sssr_davb_est(d)
% 
% Input arguments:
%     d                   array of experimental inputs and responses
%         d(:,1)              cues to subject (binary: 0 or 1)
%         d(:,2)              outcomes (binary: 0 or 1)
%         d(:,3)              reaction times (milliseconds) -- UNUSED IN THIS VERSION --
%         d(:,4)              responses (binary: 0 or 1; other values: irregular trial)
%
% Output:
%     est.u              input to subject
%         est.u.q            cue (1=high tone, 0=low tone)
%         est.u.o            outcome (1=face/high reward, 0=house/low reward)
%         est.u.x1           congruency of outcome (1=congruent, 0=incongruent)
%     est.y              observed responses
%         est.y.c            choice
%         est.y.rt           reaction time
%     est.c              configuration settings used for inversion (see function bp_configure below)
%     est.p              estimated parameters
%         est.p.mu2_0        initial value of mu2
%         est.p.sa2_0        initial value of sigma2
%         est.p.mu3_0        initial value of mu3
%         est.p.sa3_0        initial value of sigma3
%         est.p.ka           kappa
%         est.p.om           omega
%         est.p.th           theta
%         est.p.ze           zeta
%     est.la:            evolution of sufficient statistics of Gaussian approximate posteriors
%                        and other quantities of interest
%         est.la.mu1         mu1
%         est.la.mu2         mu2
%         est.la.sa2         sigma2
%         est.la.mu3         mu3
%         est.la.sa3         sigma3
%         est.la.sa1hat      precision of prediction at the 1st level
%         est.la.sa2hat      precision of prediction at the 2nd level
%         est.la.sa3hat      precision of prediction at the 3rd level
%         est.la.w2          weighting factor of informational and environmental uncertainty at the 2nd level
%         est.la.da1         prediction error at the 1st level
%         est.la.da2         prediction error at the 2nd level
%     est.posterior:    posterior as output from DAVB
%     est.out:          further output from DAVB
%
% Example:
%     est = bp_sssr_davb_est(data)
%     bp_plot_la(est)

% Christoph Mathys ETHZ/UZH
% -------------------------------------------------------------------------
function r = bp_sssr_davb_est(d)
    
    % ---------------------------------------------------------------------
    % Estimate parameters kappa, omega, theta, and zeta
    % as well as initial values of representations
    % (i.e. sufficient statistics of Gaussian approximate
    % posteriors on states x1, x2, and x3 of environment)
    % mu1, mu2, sigma2, mu3, sigma3
    
    % Store input r.u and subject responses r.y in newly
    % initialized structure r
    r = bp_data(d);
    
    % Estimate mode of posterior parameter distribution
    % and simulate corresponding representations
    tic
    r = bp_mode_est(r);
    toc
    
return;

% -------------------------------------------------------------------------
function c = bp_configure

c = struct;

% Regularization factor
c.rf = 0;

% Upper bound for kappa and theta (lower bound is always zero)
c.kaub = 6;
c.thub = 0.005;

% Sufficient statistics of Gaussian parameter priors

% Initial mu2
c.mu2_0mu = 0;
c.mu2_0sa = 1; % 1e-10

% Initial sigma2
c.logsa2_0mu = log(1);
c.logsa2_0sa = 1; % 1e-10

% Initial mu3 is gauge fixed to 1
c.mu3_0mu = 1;
c.mu3_0sa = 0;

% Initial sigma3
c.logsa3_0mu = log(1);
c.logsa3_0sa = 1; % 1e-10

% Kappa
c.logitkamu = 0;
c.logitkasa = 3^2;

% Omega
c.ommu = -4;
c.omsa = 0;

% Theta
c.logitthmu = 0;
c.logitthsa = 10^2;

% Zeta
c.logzemu = log(48);
c.logzesa = 1;

return;

% -------------------------------------------------------------------------
function r = bp_data(d)

% Initialize data structure to be returned
r = struct;

% Input: cues and outcomes
r.u.q  = d(:,1)';
r.u.o  = d(:,2)';

% Calculate congruency of input
% THIS IS THE STIMULUS (OR FUNCTION OF SEVERAL STIMULI) THE SUBJECT LEARNS FROM
r.u.x1  = (r.u.q == r.u.o);

% Responses: choices and reaction times
r.y.c  = d(:,4)';
r.y.rt = d(:,3)';

% Calculate congruency of choices
r.y.cc  = (r.u.q == r.y.c);

% Determine irregular trials
irr = [];
for k = 1:length(r.y.c)
    if r.y.c(k) ~= 0 && r.y.c(k) ~= 1
        irr = [irr, k];
    end
end

r.irr = irr;

if length(irr) == 0
    irrout = 'none';
else
    irrout = irr;
end
disp(['Irregular trials: ', num2str(irrout)])
    
return;

% -------------------------------------------------------------------------
function r = bp_mode_est(r)

% Read configuration
r.c = bp_configure;

% Get input
y = r.y.cc;
u = r.u.x1;

% Deal with irregular trials
isYout = zeros([1,length(y)]);
isYout(r.irr) = 1;

% Invert model
% i.e., OBSERVE THE OBSERVER:
f_fname = @f_VBvolatile;
g_fname = @g_VBvolatile;
dim.n = 7;
dim.n_theta = 3;
dim.n_phi = 1;
priors.muX0 = [0.5;
               r.c.mu2_0mu;
               r.c.logsa2_0mu;
               r.c.mu3_0mu;
               r.c.logsa3_0mu;
               r.c.mu2_0mu;
               r.c.logsa2_0mu];
priors.SigmaX0 = diag([0;
                    r.c.mu2_0sa;
                    r.c.logsa2_0sa;
                    r.c.mu3_0sa;
                    r.c.logsa3_0sa;
                    r.c.mu2_0sa;
                    r.c.logsa2_0sa]);
priors.muTheta = [r.c.logitkamu;
                  r.c.ommu;
                  r.c.logitthmu];
priors.SigmaTheta = diag([r.c.logitkasa;
                    r.c.omsa;
                    r.c.logitthsa]);
priors.muPhi = r.c.logzemu;
priors.SigmaPhi = r.c.logzesa;
priors.a_alpha = Inf;
priors.b_alpha = 0;
options.priors = priors;
options.binomial = 1;
options.TolFun = 1e-4;
options.inF.rf = r.c.rf;
options.inF.kaub = r.c.kaub;
options.inF.thub = r.c.thub;
options.isYout = isYout;
options.DisplayWin = 1;
options.verbose = 0;
[posterior,out] = VBA_NLStateSpaceModel(y,u,f_fname,g_fname,dim,options);

% Gather output
r.F = out.F;

r.p.mu2_0 = posterior.muX0(2);
r.p.sa2_0 = exp(posterior.muX0(3));
r.p.mu3_0 = posterior.muX0(4);
r.p.sa3_0 = exp(posterior.muX0(5));
r.p.ka      = sgm(posterior.muTheta(1),r.c.kaub);
r.p.om      = posterior.muTheta(2);
r.p.th      = sgm(posterior.muTheta(3),r.c.thub);
r.p.ze      = exp(posterior.muPhi(1));
r.p.pvec = [r.p.mu2_0, r.p.sa2_0, r.p.mu3_0, r.p.sa3_0, r.p.ka, r.p.om, r.p.th, r.p.ze];

% Store representations at mode
r.la = bp_la(u, r.p.pvec, r.c.rf);

% Store DAVB output
r.posterior = posterior;
r.out = out;

% Display estimated belief
%unwrapVBvolatileOTO(posterior,out);

% Display diagnostics
%VBA_ReDisplay(posterior, out)

% Print mode
disp(r.p)

return;
