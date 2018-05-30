function [Iphi,SigmaPhi,deltaMuPhi,suffStat] = VBA_Iphi_UNL(phi,y,posterior,suffStat,dim,u,options)
% Gauss-Newton update of the obs. parameters, for un-normalized likelihoods


if options.DisplayWin % Display progress
    try
        STR = 'VB Gauss-Newton on observation parameters... ';
        set(options.display.hm(1),'string',STR);
        set(options.display.hm(2),'string','0%');
        drawnow
    end
end

%  Look-up which evolution parameter to update
indIn = options.params2update.phi;

% Preallocate intermediate variables
Q = options.priors.SigmaPhi(indIn,indIn);
iQ = VBA_inv(Q);
muPhi0 = options.priors.muPhi;
Phi = muPhi0;
Phi(indIn) = phi;
dphi0 = muPhi0-Phi;

% inverse temperature
beta = posterior.a_sigma./posterior.b_sigma;

LL = 0;
dLLdP = zeros(1,options.dim.n_phi);
d2LLdP2 = zeros(options.dim.n_phi,options.dim.n_phi);
Ey = zeros(options.dim.p,options.dim.n_t);
Vy = zeros(options.dim.p,options.dim.n_t);
div = 0;

%--- Loop over time series ---%
for t=1:dim.n_t
    
    if ~options.isYout
        
        % evaluate re-normalized likelihood (as well as gradients and Hessians)
        [LLt,dLLdXt,dLLdPt,d2LLdX2t,d2LLdP2t,Ey(:,t),Vy(:,t)] = VBA_evalAL([],Phi,beta,u(:,t),y(:,t),options);
        
        LL = LL + LLt;
        dLLdP = dLLdP + dLLdPt;
        d2LLdP2 = d2LLdP2 + d2LLdP2t;
        
        % Accelerate divergent update
        if VBA_isWeird ({LL, dLLdP, d2LLdP2, Ey, Vy})
            div = 1;
            break
        end
        
    end
    
    % Display progress
    if mod(t,dim.n_t./10) < 1
        if  options.DisplayWin
            try
                set(options.display.hm(2),'string',[num2str(floor(100*t/dim.n_t)),'%']);
                drawnow
            end
        end
    end
    
    % Check gradients?
    if options.checkGrads
        mayPause = 0;
        if ~isempty(Phi)
            dLLdPt2 = VBA_numericDiff(@VBA_evalAL,2,[],Phi,beta,u(:,t),y(:,t),options);
            if ~ VBA_isWeird (dLLdPt2)
                [hf] = VBA_displayGrads(dLLdPt',dLLdPt2,'Gradients wrt parameters',options.g_fname,'g');
                mayPause = 1;
            else
                VBA_disp('VBA check_grads: Warning: weird numerical gradients!!!')
            end
            d2LLdP2t2 = VBA_numericDiff(@numericDiff,4,@VBA_evalAL,2,[],Phi,beta,u(:,t),y(:,t),options);
            if ~ VBA_isWeird (dLLdPt2)
                [hf] = VBA_displayGrads(d2LLdP2t,d2LLdP2t2,'Hessians wrt parameters',options.g_fname,'g');
                mayPause = 1;
            else
                VBA_disp('VBA check_grads: Warning: weird numerical gradients!!!')
            end
        end
        if mayPause
            pause
            close(setdiff(hf,0))
        end
    end
    
    
end


% Display progress
if options.DisplayWin
    try
        set(options.display.hm(2),'string','OK');
        drawnow
    end
end

% posterior covariance matrix
iSigmaPhi = iQ - d2LLdP2(indIn,indIn);
SigmaPhi = VBA_inv(iSigmaPhi);

% mode
tmp = iQ*dphi0(indIn) + dLLdP(indIn)';
deltaMuPhi = SigmaPhi*tmp;

% variational energy
Iphi = -0.5.*dphi0(indIn)'*iQ*dphi0(indIn) + LL;
if VBA_isWeird ({Iphi, SigmaPhi})
    Iphi = -Inf;
end

% update sufficient statistics
suffStat.logL = LL;
suffStat.gx = Ey; % model post-dicted observations
suffStat.dy = y - suffStat.gx;
suffStat.vy = Vy;
suffStat.dphi = dphi0;
suffStat.div = div;

