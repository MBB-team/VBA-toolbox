function [Iphi,SigmaPhi,deltaMuPhi,suffStat] = VBA_Iphi_extended(phi,y,posterior,suffStat,dim,u,options)
% Gauss-Newton update of the observation parameters
% !! When the observation function is @VBA_odeLim, this Gauss-Newton update
% actually implements a gradient ascent on the variational energy of the
% equivalent deterministic DCM.
% !! All mean-field perturbation terms are ignored (options.ignoreMF=true)


if options.DisplayWin % Display progress
    if isequal(options.g_fname,@VBA_odeLim)
        STR = 'VB Gauss-Newton on observation/evolution parameters... ';
    else
        STR = 'VB Gauss-Newton on observation parameters... ';
    end
    set(options.display.hm(1),'string',STR);
    set(options.display.hm(2),'string','0%');
    drawnow
end

%  Look-up which evolution parameter to update
indIn = options.params2update.phi;


% Preallocate intermediate variables
iQy = options.priors.iQy;
Q = options.priors.SigmaPhi(indIn,indIn);
iQ = VBA_inv(Q);
muPhi0 = options.priors.muPhi;
Phi = muPhi0;
Phi(indIn) = phi;
dphi0 = muPhi0-Phi;
dy = zeros(dim.p,dim.n_t);
vy = zeros(dim.p,dim.n_t);
gx = zeros(dim.p,dim.n_t);
ddydphi = 0;
d2gdx2 = 0;

dy2 = zeros(1,numel(options.sources));
logL = zeros(1,numel(options.sources));

if isequal(options.g_fname,@VBA_odeLim)
    clear VBA_odeLim
    muX = zeros(options.inG.old.dim.n,dim.n_t);
    SigmaX = cell(dim.n_t,1);
end
div = 0;

isTout = all (options.isYout, 1);

%--- Loop over time series ---%
for t=1:dim.n_t
    % evaluate observation function at current mode
    [gx(:,t),dG_dX,dG_dPhi] = VBA_evalFun('g',posterior.muX(:,t),Phi,u(:,t),options,dim,t);
        
    % store states dynamics if ODE mode
    if isequal(options.g_fname,@VBA_odeLim)
        % get sufficient statistics of the hidden states from unused i/o in
        % VBA_evalFun.
        muX(:,t) = dG_dX.xt;
        SigmaX{t} = dG_dX.dx'*posterior.SigmaPhi*dG_dX.dx;  
    end
    
    if ~isTout(t)
        
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
         
            [ddydphi_t,d2gdx2_t,logL_t,dy(idx_obs,t),dy2_t,vy(idx_obs,t)] = VBA_get_dL(gx(idx_obs,t),dG_dPhi(:,idx_obs),y(idx_obs,t),options.sources(si).type,iQyt,sigmaHat);
                       
            % aggregate
            ddydphi = ddydphi + ddydphi_t;
            d2gdx2  = d2gdx2 + d2gdx2_t;
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

   end 
    
    % Display progress
    if mod(t,dim.n_t./10) < 1
        if options.DisplayWin
            set(options.display.hm(2),'string',[num2str(floor(100*t/dim.n_t)),'%']);
            drawnow
        end
    end
    
    
    
    % Accelerate divergent update
    if VBA_isWeird ({dy2, dG_dPhi, dG_dX})
        div = 1;
        break
    end
    
    
end


% Display progress
if options.DisplayWin
    set(options.display.hm(2),'string','OK');
    drawnow
end

% posterior covariance matrix

iSigmaPhi = iQ + d2gdx2(indIn,indIn);
SigmaPhi = VBA_inv(iSigmaPhi);

% mode
tmp = iQ*dphi0(indIn) + ddydphi(indIn);
deltaMuPhi = SigmaPhi*tmp;

% variational energy
Iphi = -0.5.*dphi0(indIn)'*iQ*dphi0(indIn) + sum(logL);
if VBA_isWeird ({Iphi, SigmaPhi}) || div
    Iphi = -Inf;
end

% update sufficient statistics
suffStat.Iphi = Iphi;
suffStat.gx = gx;
suffStat.dy = dy;
suffStat.dy2 = dy2;
suffStat.logL = logL;
suffStat.vy = vy;
suffStat.dphi = dphi0;
if isequal(options.g_fname,@VBA_odeLim)
    suffStat.muX = muX;
    suffStat.SigmaX = SigmaX;
end
suffStat.div = div;

end



