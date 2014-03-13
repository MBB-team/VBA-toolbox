% Fits the Belief-Precision model to observed data. Top-level function.
% 
% Usage:
%     est = bp_ptrt_davb_est('data.txt', 'respmod')
% 
% Input:
%     'data.txt'             Name of a data text file organized in four columns
%         1st column         Trial number
%         2nd column         Cue (1 or 2)
%         3rd column         Target (1 or 2)
%         4th column         Reaction time in milliseconds (0 = error)
%     'respmod'              The type of response model to be used. Possible values
%                            are 'precision', 'surprise', or 'belief'
%
% Output:
%     est.u              input to subject
%         est.u.q            cue (1 or 2)
%         est.u.t            target (1 or 2)
%         est.u.x1           validity (i.e. congruency of cue and target) (1=valid, 0=invalid)
%     est.y              observed responses
%         est.y.rt           reaction time (0 = error)
%     est.irr            trial numbers of irregular trials
%     est.c              configuration settings used for inversion (see function bp_configure below)
%     est.p              estimated parameters
%         est.p.mu2_0        initial value of mu2
%         est.p.sa2_0        initial value of sigma2
%         est.p.mu3_0        initial value of mu3
%         est.p.sa3_0        initial value of sigma3
%         est.p.ka           kappa
%         est.p.om           omega
%         est.p.th           theta
%         est.p.ze1          zeta1
%         est.p.ze2          zeta2
%         est.p.ze3          zeta3
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
%     est = bp_ptrt_davb_est('data.txt', 'precision');
%     bp_plot_la(est)

% Christoph Mathys ETHZ/UZH
% -------------------------------------------------------------------------
function r = bp_ptrt_davb_est(datafile, respmod)
    
    % ---------------------------------------------------------------------
    % Estimate parameters omega, kappa, theta, zeta1, zeta2,
    % and zeta3, as well as initial values of representations
    % (i.e. sufficient statistics of Gaussian approximate
    % posteriors on states x1, x2, and x3 of environment)
    % mu1, mu2, sigma2, sigma3
    
    % Store input r.u and subject responses r.y in newly
    % initialized structure r
    d = dlmread(datafile);
    r = bp_data(d);
    
    % Estimate mode of posterior parameter distribution
    % and simulate corresponding representations
    tic
    r = bp_mode_est(r, respmod);
    toc
    
return;

% -------------------------------------------------------------------------
function c = bp_configure(respmod)

c = struct;

% Response model
c.respmod = respmod;

% Regularization factor
c.rf = 0;

% Upper bound for kappa and theta (lower bound is always zero)
c.kaub = 1.4;
c.thub = 1;

% Sufficient statistics of Gaussian parameter priors

% Initial mu2
c.mu2_0mu = 0.85; % corresponds to sgm(mu2) = 0.70
c.mu2_0sa = 1e-10; %(voher 1e-10)

% Initial sigma2
c.logsa2_0mu = log(0.06);
c.logsa2_0sa = 1e-10;

% Initial mu3 is gauge fixed to 1
c.mu3_0mu = 1;
c.mu3_0sa = 0;

% Initial sigma3
c.logsa3_0mu = log(4);
c.logsa3_0sa = 1e-10;

% Kappa
c.logitkamu = 0;
c.logitkasa = 10^2;

% Omega
c.ommu = -4;
c.omsa = 10;

% Theta
c.logitthmu = 0;
c.logitthsa = 10^2;

% Zeta1
c.logze1mu = log(0.0052);
c.logze1sa = 1e-1;

% Zeta2
c.logze2mu = log(0.0006); % Corresponds to sigmoid parameter = 1
c.logze2sa = 1e-10;

% Zeta3
c.ze3mu = 0.001; % Corresponds to inverse RT 2noise of 0.001 1/ms
c.ze3al = 10;

return;

% -------------------------------------------------------------------------
function r = bp_data(d)

% Initialize data structure to be returned
r = struct;

% Input: cues and outcomes
r.u.q  = d(:,2)';
r.u.t  = d(:,3)';

% Calculate congruency of input
% THIS IS THE STIMULUS (OR FUNCTION OF SEVERAL STIMULI) THE SUBJECT LEARNS FROM
r.u.x1  = (r.u.q == r.u.t);

% Responses: reaction times
r.y.rt = d(:,4)';

% Determine irregular trials
irr =  d(d(:,4)==0,1)';
r.irr = irr;

if length(irr) == 0
    irrout = 'none';
else
    irrout = irr;
end
disp(['Irregular trials: ', num2str(irrout)])

return;

% -------------------------------------------------------------------------
function r = bp_mode_est(r, respmod)

% Read configuration
r.c = bp_configure(respmod);

% Get input
y = 1./r.y.rt; % INVERSE reaction times!
y(r.irr) = 0;
u = r.u.x1;

% Deal with irregular trials
iQy = cell(length(u),1);
iQy(:) = {1};
iQy(r.irr) = {0};

% Invert model
% i.e., OBSERVE THE OBSERVER:
f_fname = @f_VBvolatile;
g_fname = @g_VBvolatile;
dim.n = 7;
dim.n_theta = 3;
dim.n_phi = 2;
priors.iQy = iQy;
priors.muX0 = [0.5; %sgm(r.c.mu2_0mu,1);
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
priors.muPhi = [r.c.logze1mu;
                r.c.logze2mu];
priors.SigmaPhi = diag([r.c.logze1sa;
                    r.c.logze2sa]);
priors.a_alpha = Inf;
priors.b_alpha = 0;
priors.a_sigma = r.c.ze3al;
priors.b_sigma = r.c.ze3mu/r.c.ze3al;
options.priors = priors;
options.binomial = 0;
options.TolFun = 1e-4;
options.inF.rf = r.c.rf;
options.inF.kaub = r.c.kaub;
options.inF.thub = r.c.thub;
options.inG.respmod = respmod;
options.DisplayWin = 1;
options.verbose = 1;
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
r.p.ze1     = exp(posterior.muPhi(1));
r.p.ze2     = exp(posterior.muPhi(2));
r.p.ze3     = (posterior.a_sigma-1)*posterior.b_sigma;
r.p.pvec = [r.p.mu2_0, r.p.sa2_0, r.p.mu3_0, r.p.sa3_0, r.p.ka, r.p.om, r.p.th, r.p.ze1, r.p.ze2, r.p.ze3];

% Store representations at mode
r.la = bp_la(r.u.x1, r.p.pvec, r.c.rf);

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
