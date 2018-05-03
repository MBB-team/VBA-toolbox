function [suffStat,posterior] = VBA_check_errors(y,u,options)
% dummy diagnostic of model specification

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
dy2 = zeros(1,numel(options.sources));
logL = zeros(1,numel(options.sources));

if  sum([options.sources(:).type]==0)>0
    iQy = options.priors.iQy; 
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
        if VBA_isWeird (gx(find (options.isYout(:, t) == 0), t)) 
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


   % accumulate gradients, hessian and likelyhood
    gsi = find([options.sources(:).type]==0) ;
    for si = 1:numel(options.sources)
        % compute source contribution
        idx_obs_all = options.sources(si).out;
        is_obs_out = find(options.isYout(idx_obs_all,t)==0);
        idx_obs = idx_obs_all(is_obs_out);
       
        if ~isempty(idx_obs) 
            
            sigmaHat=0;
            iQyt=[];
            if options.sources(si).type==0
                gi = find(si==gsi) ;
                sigmaHat = posterior.a_sigma(gi)./posterior.b_sigma(gi);
                iQyt = iQy{t,gi};
                iQyt=iQyt(is_obs_out,is_obs_out);
            end
            
            [~,~,logL_t,dy(idx_obs,t),dy2_t,vy(idx_obs,t)] = VBA_get_dL(gx(idx_obs,t),dG_dPhi(:,idx_obs),y(idx_obs,t),options.sources(si).type,iQyt,sigmaHat);
            
            % aggregate
            dy2(si) = dy2(si) + dy2_t;
            logL(si) = logL(si) + logL_t;
            
             if options.sources(si).type==0
                V = dG_dPhi(:,idx_obs)'*posterior.SigmaPhi*dG_dPhi(:,idx_obs) ;
                if dim.n > 0
                    V = V + dG_dX(:,idx_obs)'*posterior.SigmaX.current{t}*dG_dX(:,idx_obs);
                end
                vy(idx_obs,t) = vy(idx_obs,t) + diag(V);
            end
                 
        end
        
        
        
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
    if VBA_isWeird ({dy(find(options.isYout(:,t)==0),t), dG_dPhi(:,find(options.isYout(:,t)==0)), dG_dX})
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
suffStat.dy2 = dy2;
suffStat.logL = logL;





