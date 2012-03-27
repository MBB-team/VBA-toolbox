% % demo simulating hrf with SPM functions


clear variables
n = 1;
 
T      = 64;   % seconds
TR     = 2;    % seconds
U.dt   = TR/32; % micro-time seconds
U.name = 'dummy';
U.u    = zeros(T/U.dt,1);
t      = (1:T/U.dt)*U.dt;
i      = find(t) > 8 & t < 32;
U.u(i) = -1;
 
a = eye(n);
b = zeros(n);
c = zeros(n,1);
c(1) = 1;
d = zeros(n,n,n);
[pE,pC] = spm_dcm_fmri_priors(a,b,c,d);
pE.C(1) = 1;
x = sparse(n,5);
 
 
M.states = 1:n;
M.nsteps = round(TR);
M.delays = TR*ones(n,1)/2;
M.TE = 0.04;
M.f  = 'spm_fx_fmri';
M.g  = 'spm_gx_fmri';
M.x  = x;
M.pE = pE;
M.pC = pC;
M.m  = size(U.u,2);
M.n  = size(x(:),1);
M.l  = size(x,1);
M.ns = T/TR;
M.W  = 1e6;
 
SPM_INT = {
    'spm_int', ...
    'spm_int_B', ...
    'spm_int_D', ...
    'spm_int_J', ...
    'spm_int_L', ...
    'spm_int_ode', ...
    'spm_int_sde'};
 
ha = axes('parent',figure,'nextplot','add');
col = { 'r', 'g' , 'b' , 'c' , 'm' , 'k' , 'y' };
for i=1:length(SPM_INT)
    y{i} = feval(SPM_INT{i},M.pE,M,U);
    t    = U.dt*length(U.u)*[1:length(y{i})]/(length(y{i}));
    plot(ha,t,y{i},col{i});
end
legend(SPM_INT,'interpreter','none')
grid(ha,'on')


% clear variables
% n  = 1; % number of regions
% v = 256; % number of scans
% 
% U.dt = 32/v;
% U.name = 'dummy';
% U.u  = zeros(v,1);
% U.u(2:20) = 1./U.dt;
% % U.u(2:30) = -1./U.dt; % increase inhibitory block
% 
% a = eye(n);
% b = zeros(n);
% c = zeros(n,1);
% c(1) = 1;
% d = zeros(n,n,n);
% [pE,pC] = spm_dcm_fmri_priors(a,b,c,d);
% pE.C(1) = 1;
% x = sparse(n,5);
% 
% M.TE = 0.04;
% % M.delays = 0*ones(n,1);  % NB: including this provokes error in spm_int!
% M.states = 1:n;
% M.f = 'spm_fx_fmri';
% M.g = 'spm_gx_fmri';
% M.x = x;
% M.pE = pE;
% M.pC = pC;
% M.m = size(U.u,2);
% M.n = size(x(:),1);
% M.l = size(x,1);
% M.ns = v;
% M.W  = speye(n*5,n*5)*1e16;
% 
% SPM_INT = {
%     'spm_int',...
%     'spm_int_B',...
%     'spm_int_D',...
%     'spm_int_df',...
%     'spm_int_E',...
%     'spm_int_J',...
%     'spm_int_L',...
%     'spm_int_LP',...
%     'spm_int_ode',...
%     'spm_int_sde'};
% 
% ha = axes('parent',figure,'nextplot','add');
% warning off
% out = [];
% for i=1:length(SPM_INT)
%     try
%         [y(:,i)] = feval(SPM_INT{i},M.pE,M,U);
%         if isweird(y(:,i))
%             out = [out;i];
%             disp([SPM_INT{i},' gave NaN or Inf time series.'])
%         end
%     catch
%         disp([SPM_INT{i},' did not work (missing argument?).'])
%         out = [out;i];
%     end
% end
% SPM_INT(out) = [];
% y(:,out) = [];
% 
% plot(ha,y)
% legend(SPM_INT,'interpreter','none')
% grid(ha,'on')
% 
% 
% 
