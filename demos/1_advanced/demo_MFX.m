function demo_MFX
% this demo exemplifies the use of mixed-effects analysis in VBA

k = 1;

ns = 50; % #subjects
dim.n_phi = 1;
dim.n = 0;
dim.n_theta = 0;
dim.p = 1; % data dim (within-subject)
dim.n_t = 10; % trials (within-subject)

% simulate MFX


phi_muGroup = ones(1, k);
phi_SigmaGroup = eye(k);


y = cell(ns,1);
sigmas = 0.1 * ones(ns,1); % within-subject residual variance
g_fname = @g_mfx;
f_fname = [];
for i=1:ns
    % draw within-subject effects from population distribution
    phi(:,i) = VBA_random ('Gaussian', phi_muGroup, phi_SigmaGroup);
    u{i} = randn(k,dim.n_t);
    options{i}.DisplayWin = 0;
    options{i}.verbose = 0;
    options{i}.dim = dim;
    [y{i}] = VBA_simulate (dim.n_t,f_fname,g_fname,[],phi(:,i),u{i},Inf,sigmas(i),options{i});
end

priors_group = struct();
priors_group.SigmaPhi = 10 * eye(dim.n_phi);


% % TEST CASES (to comment/uncomment)
% % 1. fix population mean to 0 for phi(1)
% priors_group.SigmaPhi(1,1) = 0;
% % 2. fixed-effect over the population for phi(1)
% priors_group.a_vPhi = Inf;
% priors_group.b_vPhi = 0;

[p_sub, o_sub, p_group, o_group] = VBA_MFX(y,u,f_fname,g_fname,dim,struct('priors',priors_group));

% extract within-subject parameter estimates
% without MFX-type priors
posterior_muPhi_RFX = [o_group.initVBA.p_sub.muPhi];
% with MFX-type priors
posterior_muPhi_MFX = [p_sub.muPhi];
posterior_muPhi_MFX_group = p_group.muPhi;

% plot
VBA_figure();
boxplot([posterior_muPhi_RFX', posterior_muPhi_MFX'], 'Whisker',Inf);
box off
hold on;
plot([.5 3.5] , mean(phi) * [1 1],'g');
xlabel('method');
set(gca,'XTickLabel',{'RFX','MFX'})
ylabel('parameter estimate')
legend('true group average', 'Location','north')


end

%% ########################################################################

function [fx] = f_mfx(x,P,u,in)
fx =   VBA_sigmoid(P(2) * P(1)) * x + P(2) * u - P(1);
end

function  [gx, dgdx, dgdp] = g_mfx(x,P,u,in)
gx = P'*u;
dgdx = [];
dgdp = u;
end
