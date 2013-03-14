function [suffStat,posterior] = VBA_check_errors_extended(y,u,options)
% Function merging VBA_getsuffstat and the content of VBA_Initialize
% regarding deterministic DCM

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

if ~options.binomial || sum([options.sources(:).type]==0)>1
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
            VB_disp('Output dimensions of either the observation or evolution function are incorrect',options)
        end
        VBA_disp(' ',options)
        posterior = [];
        return
    end
    
 
    
%     %- check gaussian sources
%     gsi = find([options.sources(:).type]==0) ;
%     for i=1:length(gsi)   
%         si = gsi(i) ; % source number
%         s_out = options.sources(si).out ; % corresponding observation indices
%         % Get precision parameters
%         sigmaHat = posterior.a_sigma(si)./posterior.b_sigma(si);
%          % error terms
%         s_yin = s_out(~options.isYout(s_out,t));
% %         if length(s_yin) == length(s_out)
%             dy(s_yin,t) = y(s_yin,t) - gx(s_yin,t);
%             dy2_t = dy(s_yin,t)'*iQy{t,si}*dy(s_yin,t);
%             dy2(gsi(si)) = dy2(gsi(si)) + dy2_t;
%             logL(si) = 0;% logL(si) -.5*sigmaHat*dy2_t;
% %         end
%         % Predictive density (data space)
%         V = dG_dPhi(:,s_out)'*posterior.SigmaPhi*dG_dPhi(:,s_out) + (1./sigmaHat).*VB_inv(iQy{t,si},[]);
%         vy(s_out,t) = diag(V);
%     end
%     %- check binomial sources
%     bsi = find([options.sources(:).type]==1) ;
%     for i=1:length(bsi)
%         si = bsi(i) ; % source number
%         s_out = options.sources(si).out ; % corresponding observation indices
%         % fix numerical instabilities
%         gx(s_out,t) = checkGX_binomial(gx(s_out,t)); 
%         % predicted variance over binomial data
%         vy(s_out,t) = gx(s_out,t).*(1-gx(s_out,t));
%         s_yin = s_out(~options.isYout(s_out,t));
%         % accumulate log-likelihood
%         logL(si) = logL(si) + y(s_yin,t)'*log(gx(s_yin,t)) + (1-y(s_yin,t))'*log(1-gx(s_yin,t));
%         % prediction error
%         dy(s_yin,t) = y(s_yin,t) - gx(s_yin,t);
%     end
%     %- check multinomial sources
%     msi = find([options.sources(:).type]==2) ;
%     for i=1:length(msi)
%         si = msi(i) ; % source number
%         s_out = options.sources(si).out ; % corresponding observation indices
%         % fix numerical instabilities
%         gx(s_out,t) = checkGX_binomial(gx(s_out,t));
%         % predicted variance over binomial data
%         vy(s_out,t) = gx(s_out,t).*(1-gx(s_out,t));
%         s_yin = s_out(~options.isYout(s_out,t));
%         % accumulate log-likelihood
%         logL(si) = logL(si) + log(gx(s_yin,t))'*y(s_yin,t) ;
%         % prediction error
%         dy(s_yin,t) = y(s_yin,t) - gx(s_yin,t);
%         
%     end
%%
   % accumulate gradients, hessian and likelyhood
    gsi = find([options.sources(:).type]==0) ;
    for si = 1:numel(options.sources)
        % compute source contribution
        idx_obs_all = options.sources(si).out;
        idx_obs = idx_obs_all(options.isYout(idx_obs_all,t)==0);
       
        if ~isempty(idx_obs) 
            
            sigmaHat=0;
            iQyt=[];
            if options.sources(si).type==0
                gi = find(si==gsi) ;
                sigmaHat = posterior.a_sigma(gi)./posterior.b_sigma(gi);
                iQyt = iQy{t,gi};
            end
            
            [~,~,logL_t,dy(idx_obs,t),dy2_t,vy(idx_obs,t)] = VBA_get_dL(gx(idx_obs,t),dG_dPhi(:,idx_obs),y(idx_obs,t),options.sources(si).type,iQyt,sigmaHat);
                       
            % aggregate
            dy2(si) = dy2(si) + dy2_t;
            logL(si) = logL(si) + logL_t;

        end
        
        
        
    end

    % include parameter variance in predictive density of data
    Vsi = [options.sources(gsi).out];
    if ~isempty(Vsi)
    	V = dG_dPhi(:,Vsi)'*posterior.SigmaPhi*dG_dPhi(:,Vsi) ;
        if dim.n > 0
            V = V + dG_dX(:,Vsi)'*posterior.SigmaX.current{t}*dG_dX(:,Vsi);
        end
        vy(Vsi,t) = vy(Vsi,t) + diag(V);
    end

 %%
    
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


suffStat.dy2 = dy2;
suffStat.logL = logL;





