function [suffStat,posterior] = VBA_check_errors(y,u,options)
% dummy diagnostic of model specification

if options.extended
    [suffStat,posterior] = VBA_check_errors_extended(y,u,options);
    return
end

[suffStat] = VBA_getSuffStat(options,[],0);
dim = options.dim;
posterior = options.priors;
posterior.muX = sparse(0,dim.n_t);

% watch out about display setting
options.DisplayWin = 0;
if ~options.OnLine && options.verbose
    fprintf(1,'Deriving prior''s sufficient statistics ...')
    fprintf(1,'%6.2f %%',0)
end


% Preallocate intermediate variables
Phi = options.priors.muPhi;
dy = zeros(dim.p,dim.n_t);
vy = zeros(dim.p,dim.n_t);
gx = zeros(dim.p,dim.n_t);
if ~options.binomial
    iQy = options.priors.iQy;
    dy2 = 0;
    % Get precision parameters
    sigmaHat = posterior.a_sigma./posterior.b_sigma;
else
    logL = 0;
end
if isequal(options.g_fname,@VBA_odeLim)
    clear VBA_odeLim
    muX = zeros(options.inG.old.dim.n,dim.n_t);
    SigmaX = cell(dim.n_t,1);
end
div = 0;

%--- Loop over time series ---%
for t=1:dim.n_t
    
    % evaluate observation function at current mode
    try
        [gx(:,t),dG_dX,dG_dPhi] = VBA_evalFun('g',posterior.muX(:,t),Phi,u(:,t),options,dim,t);
        if isweird(gx(:,t))
            VBA_disp('',options)
            VBA_disp('Error: could not initialize VB scheme: model generates NaN or Inf!',options)
            posterior = [];
            return
        end
    catch ME
        name = ME.stack.name; % name of the function in which error occured
        line = ME.stack.line; % number of the line in which error occured
        file = ME.stack.file; % location of the function in which the error occured
        link = ['<a href = "matlab: open ',name,'">',name,'</a>'];
        VBA_disp(' ',options)
        VBA_disp(['Error: could not initialize VB scheme : check function "',link,'" at line ', num2str(line)],options)
        VBA_disp('---------------',options)
        VBA_disp('CAUSE :',options)
        VBA_disp([ME.message],options);
        fid = fopen([name,'.m'],options);
        for l = 1:line
            codeline = fgets(fid);
        end
        fclose(fid);
        VBA_disp(codeline,options)
        VBA_disp('---------------',options)
        if isequal(ME.message,'Subscripted assignment dimension mismatch.')
            VBA_disp('Output dimensions of either the observation or evolution function are incorrect',options)
        end
        VBA_disp(' ',options)
        posterior = [];
        return
    end
    
    % remove irregular trials
    yin = find(~options.isYout(:,t));
    if ~options.binomial
        % error terms
        dy(:,t) = y(:,t) - gx(:,t);
        dy2 = dy2 + dy(:,t)'*iQy{t}*dy(:,t);
        % Predictive density (data space)
        V = dG_dPhi'*posterior.SigmaPhi*dG_dPhi + (1./sigmaHat).*VBA_inv(iQy{t},[]);
        vy(:,t) = diag(V);
    else
        % fix numerical instabilities
        gx(:,t) = checkGX_binomial(gx(:,t));
        % predicted variance over binomial data
        vy(:,t) = gx(:,t).*(1-gx(:,t));
        % accumulate log-likelihood
        logL = logL + y(yin,t)'*log(gx(yin,t)) + (1-y(yin,t))'*log(1-gx(yin,t));
        % prediction error
        dy(yin,t) = y(yin,t) - gx(yin,t);
    end
    
    % store states dynamics if ODE mode
    if isequal(options.g_fname,@VBA_odeLim)
        % get sufficient statistics of the hidden states from unused i/o in
        % VBA_evalFun.
        muX(:,t) = dG_dX.xt;
        SigmaX{t} = dG_dX.dx'*posterior.SigmaPhi*dG_dX.dx;
    end
    
    % Display progress
    if mod(t,dim.n_t./10) < 1
        if ~options.OnLine && options.verbose
            fprintf(1,repmat('\b',1,8))
            fprintf(1,'%6.2f %%',100*t/dim.n_t)
        end
    end
    
    % Accelerate divergent update
    if isweird({dy,dG_dPhi,dG_dX})
        div = 1;
        VBA_disp(' ',options)
        VBA_disp('Error: Output of evolution or observation function is weird (nan or inf)   ...',options)
        VBA_disp(' ',options)
        posterior = [];
        break
    end
    
end

% Display progress
if ~options.OnLine  && options.verbose
    fprintf(1,repmat('\b',1,8))
    fprintf(' OK.')
    fprintf('\n')
end

% update sufficient statistics
suffStat.div = div;
suffStat.gx = gx;
suffStat.dy = dy;
suffStat.vy = vy;
if isequal(options.g_fname,@VBA_odeLim)
    suffStat.muX = muX;
    suffStat.SigmaX = SigmaX;
end
if ~options.binomial
    suffStat.dy2 = dy2;
else
    suffStat.logL = logL;
end