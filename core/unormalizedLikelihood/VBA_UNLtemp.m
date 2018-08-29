function [posterior,suffStat] = VBA_UNLtemp(y,posterior,suffStat,dim,u,options)
% VB update of the inverse tempterature of un-normalized likelihoods


if options.DisplayWin % Display progress
    try
        STR = 'VB update of the temperature hyperparameter... ';
        set(options.display.hm(1),'string',STR);
        set(options.display.hm(2),'string','0%');
        drawnow
    end
end

% Preallocate intermediate variables
a0 = posterior.a_sigma;
b0 = posterior.b_sigma;
beta0 = a0./b0;

c = 0;
d = 0;
s1 = 0;
s2 = 0;
div = 0;
Vy = zeros(options.dim.p,options.dim.n_t);

%--- Loop over time series ---%
for t=1:dim.n_t
    
    if ~options.isYout
        
        % evaluate re-normalized likelihood (as well as gradients and Hessians)
        [EUi,Ui,Zi,d2Uidx2,d2UidP2,Vy(:,t)] = VBA_evalAL2([],posterior.muPhi,beta0,u(:,t),y(:,t),options);
        
        c = c - beta0*EUi;
        d = d - Ui - 0.5*trace(posterior.SigmaPhi*d2UidP2);
        s1 = s1 + Ui;
        s2 = s2 + Zi;
        
        % Accelerate divergent update
        if VBA_isWeird({EUi,Ui,Zi,d2Uidx2,d2UidP2})
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
    
    
end


% Display progress
if options.DisplayWin
    try
        set(options.display.hm(2),'string','OK');
        drawnow
    end
end

% figure('name',['VB update of beta: beta0=',num2str(beta0)]),plot(Vy)

% posterior covariance matrix
posterior.a_sigma = options.priors.a_sigma + c;
posterior.b_sigma = options.priors.b_sigma + d;



% update sufficient statistics
b = posterior.a_sigma/posterior.b_sigma;
suffStat.logL = b*s1 - s2;
suffStat.vy = Vy;
suffStat.div = div;



