function [posterior,out] = VBA_onlineWrapper(y,u,f_fname,g_fname,dim,options)
% online VB inversion of nonlinear stochastic DCMs
% function [posterior,out] =
% VBA_onlineWrapper(y,u,f_fname,g_fname,dim,options)
% This function inverts a nonlinear state-space model 'on line', i.e.
% estimate hidden states, observation/evolution parameters and
% hyperparameters in 'real-time', by adding the time samples one after the
% other. It returns a time-dependent posterior density over all unknown
% variables. See VBA_NLStateSpaceModel.m for I/O arguments.
% NB: this online wrapper does not deal with ODE state-space models.


tStart = tic;

%------------------ Check input consistency ---------------%
% if isinf(options.priors.a_alpha) && isequal(options.priors.b_alpha,0)
%     error('Cannot invert on-line ''ODE''-like state-space models!')
%     posterior = [];
%     out = [];
%     return
% end
options.OnLine = 1;
[options,u,dim] = VBA_check(y,u,f_fname,g_fname,dim,options);

if numel(options.sources) > 1
    error('*** online inversion is not yet compatible with multisource observations');
end

if options.sources.type > 1
    error('*** online inversion is not yet compatible with multinomial observations');
end

%------------------- Initialize variables ------------------%
[suffStat] = VBA_getSuffStat(options,[],1);
suffStat.F = 0;
posterior.muX = zeros(dim.n,dim.n_t);
posterior.SigmaX.current = cell(dim.n_t,1);
if dim.n_phi >= 1
    posterior.muPhi = zeros(dim.n_phi,dim.n_t);
    posterior.SigmaPhi = cell(dim.n_t,1);
end
if dim.n_theta >= 1
    posterior.muTheta = zeros(dim.n_theta,dim.n_t);
    posterior.SigmaTheta = cell(dim.n_t,1);
end
posterior.a_alpha = zeros(1,dim.n_t);
posterior.b_alpha = zeros(1,dim.n_t);
if options.sources.type == 0
    posterior.a_sigma = zeros(1,dim.n_t);
    posterior.b_sigma = zeros(1,dim.n_t);
end

hfp = findobj('tag','VBNLSS');
if ~isempty(hfp)
    clf(hfp)
end
[options] = VBA_initDisplay(options);
if options.DisplayWin
    pos0 = get(0,'screenSize');
    set(options.display.hfp,'name','On-line version of the VB inversion',...
        'position',[0.05*pos0(3),0.05*pos0(4),0.45*pos0(3),0.9*pos0(4)])
    set(options.display.ho,'string',...
        ['Processing time point # ',num2str(1),...
        '/',num2str(dim.n_t),' ...'])
    set(options.display.hfp,'tag','online_VBNLSS')
    try,xlabel(options.display.ha(6),'time','fontsize',8);end
    try,xlabel(options.display.ha(8),'time','fontsize',8);end
    try,delete(options.display.hm(1:2));end
    OL_options = rmfield(options,'display');
else
    OL_options = options;
end
OL_dim = dim;


%-------- First data sample --------%
disp(' ')
disp('On-line VB stochastic system inversion ...')
disp(' ')
if ~options.DisplayWin
    STR_time = ['Processing time point # 1/',num2str(dim.n_t),' ...'];
    fprintf(1,STR_time)
    STR_time = [STR_time,'     '];
end

OL_options.updateX0 = 1;
OL_options.GnFigs = 0;
OL_options.DisplayWin = 0;
OL_options.isYout = options.isYout(:,1);
OL_dim.n_t=1;
OL_y = y(:,1);
OL_u = u(:,1);
[OL_posterior,OL_out] = VBA_NLStateSpaceModel(...
    OL_y,OL_u,f_fname,g_fname,OL_dim,OL_options);

[posterior,suffStat,options] = updatePSS(...
        OL_posterior,OL_out,posterior,suffStat,1,options);

posterior.muX0 = OL_posterior.muX0;
posterior.SigmaX0 = OL_posterior.SigmaX0;
suffStat.dx0 = OL_out.suffStat.dx0;
suffStat.SX0 = OL_out.suffStat.SX0;

% Display first data sample inference results
VBA_updateDisplay(posterior,suffStat,options,y(:,1),1,'precisions')
if dim.n_phi > 0
    VBA_updateDisplay(posterior,suffStat,options,y(:,1),1,'phi')
    if options.DisplayWin
        xlabel(options.display.ha(5),'time','fontsize',8)
    end
end
if dim.n > 0
    VBA_updateDisplay(posterior,suffStat,options,y(:,1),1,'X')
end
if dim.n_theta > 0
    VBA_updateDisplay(posterior,suffStat,options,y(:,1),1,'theta')
    if options.DisplayWin
        xlabel(options.display.ha(7),'time','fontsize',8)
    end
end

%-------- Loop over data samples --------%
for t =2:dim.n_t
    if options.DisplayWin
        set(options.display.ho,'string',...
            ['Processing time point # ',num2str(t),...
            '/',num2str(dim.n_t),' ...'])
    else
        back = repmat('\b',1,length(STR_time));
        fprintf(1,back)
        STR_time = ['Processing time point # ',num2str(t),...
            '/',num2str(dim.n_t),' ...'];
        fprintf(1,STR_time)
    end
    % Get last output structure
    in.out = OL_out;
    in.out.options.GnFigs = 0;
    in.out.options.DisplayWin = 0;
    % Initialize next posterior with past posterior
    in.posterior = OL_posterior;
    % Define priors for parameters with past posterior
    in.out.options.priors = OL_posterior;
    in.out.options.priors.iQy = options.priors.iQy;
    in.out.options.priors.iQx = options.priors.iQx;
    % Define initial conditions as past posterior on hidden states
    in.posterior.muX0 = OL_posterior.muX;
    in.posterior.SigmaX0 = OL_posterior.SigmaX.current{1};
    in.out.options.priors.muX0 = in.posterior.muX0;
    in.out.options.priors.SigmaX0 = in.posterior.SigmaX0;
    in.out.options.updateX0 = 0;
    % Extract current data and input
    OL_y = y(:,t);
    OL_u = u(:,t);
    % Call inversion routine on one data point
    [OL_posterior,OL_out] = VBA_NLStateSpaceModel(...
        OL_y,OL_u,[],[],[],[],in);
    % Store hidden states current posterior
    [posterior,suffStat,options] = updatePSS(...
        OL_posterior,OL_out,posterior,suffStat,t,options);
    
    % update display
    options_OL = options;
    options_OL.isYout = options_OL.isYout(:,1:t);
    VBA_updateDisplay(posterior,suffStat,options_OL,y(:,1:t),t,'precisions')
    if dim.n_phi > 0
        VBA_updateDisplay(posterior,suffStat,options_OL,y(:,1:t),t,'phi')
    end
    if dim.n > 0
        VBA_updateDisplay(posterior,suffStat,options_OL,y(:,1:t),t,'X')
    end
    if dim.n_theta > 0
        VBA_updateDisplay(posterior,suffStat,options_OL,y(:,1:t),t,'theta')
    end
end

dt = toc(tStart);
mins = floor(dt./60);
if mins >= 1
    str = ['On line inversion: complete (took ~',num2str(mins),' min).'];
else
    str = 'On line inversion: complete (took <1 min).';
end
if options.DisplayWin
    set(options.display.ho,'string',str)
else
    back = repmat('\b',1,length(STR_time));
    fprintf(1,back)
    disp(str)
end

out.F = suffStat.F(end);
out.options = options;
out.u = VBA_getU(u,options,dim,'back2micro');
out.y = y;
out.dim = dim;
out.it = dim.n_t;
out.suffStat = suffStat;
out.date = clock;
out.dt = dt;

try % recover prior predictive density
    [out.options.priors.muX,out.options.priors.SigmaX.current] = ...
        VBA_EKF(y,out.u,options.priors,dim,options,2);
end

if options.DisplayWin
    % display diagnostics
    figure(options.display.hfp)
    VBA_ReDisplay(posterior,out);
end


function [posterior,suffStat,options] = updatePSS(...
    OL_posterior,OL_out,posterior,suffStat,t,options)
posterior.muX(:,t) = OL_posterior.muX;
try
    posterior.SigmaX.current{t} = OL_posterior.SigmaX.current{1};
catch
    posterior.SigmaX.current{t} = [];
end
suffStat.SX = suffStat.SX + OL_out.suffStat.SX;
if OL_out.dim.n_phi >= 1
    posterior.muPhi(:,t) = OL_posterior.muPhi;
    posterior.SigmaPhi{t} = OL_posterior.SigmaPhi;
    suffStat.dphi(:,t) = options.priors.muPhi - posterior.muPhi(:,t);
    suffStat.Sphi = OL_out.suffStat.Sphi;
else
    posterior.muPhi = [];
    posterior.SigmaPhi = [];
    suffStat.Sphi = [];
end
if OL_out.dim.n_theta >= 1
    posterior.muTheta(:,t) = OL_posterior.muTheta;
    posterior.SigmaTheta{t} = OL_posterior.SigmaTheta;
    suffStat.dtheta(:,t) = options.priors.muTheta - posterior.muTheta(:,t);
    suffStat.Stheta = OL_out.suffStat.Stheta;
else
    posterior.muTheta = [];
    posterior.SigmaTheta = [];
    suffStat.Stheta = [];
end
try
    posterior.a_alpha(t) = OL_posterior.a_alpha;
    posterior.b_alpha(t) = OL_posterior.b_alpha;
catch
    posterior.a_alpha(t) = [];
    posterior.b_alpha(t) = [];
end
if options.sources.type == 0
    try
        posterior.a_sigma(t) = OL_posterior.a_sigma;
        posterior.b_sigma(t) = OL_posterior.b_sigma;
    catch
        posterior.a_sigma(t) = [];
        posterior.b_sigma(t) = [];
    end
end
% Store sufficient statistics
suffStat.gx(:,t) = OL_out.suffStat.gx;
suffStat.vy(:,t) = OL_out.suffStat.vy;
suffStat.dx(:,t) = OL_out.suffStat.dx;
suffStat.dy(:,t) = OL_out.suffStat.dy;
suffStat.Salpha = OL_out.suffStat.Salpha;
if options.sources.type == 0
    suffStat.Ssigma = OL_out.suffStat.Ssigma;
else
    suffStat.logL = OL_out.suffStat.logL;
end
% store integral of free energies
suffStat.F = [suffStat.F,suffStat.F(end) + OL_out.F];

% store prior predictive density
try
    options.priors.muX(:,t) = OL_out.options.priors.muX;
    options.priors.SigmaX.current{t} = OL_out.options.priors.SigmaX.current{1};
    options.priors.SigmaX.inter{t} = [];
end


