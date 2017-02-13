function [fx,indlev] =RecToMfunction(x,P,u,inF)
% recursive evolution function for k-ToM's learning rule (2x2x2 game)
% [fx,indlev] =RecToMfunction(x,P,u,inF)
% Marie Devaine wrote this function in November 2015 (comments: JD).
% A k-ToM agent with k>0 learns the unknown parameters (theta) of her
% opponent's learning and decision rules. Typically, theta includes her
% opponent's prior volatility about herself (which controls the learning
% rate), the behavioural temperature and the bias. In addition, when k>1,
% k-ToM also learns her opponent's sophistication level (k').
% If the agent knew k', then the belief update of theta would simply be
% proportional to the prediction error PE = o-P(o=1|k',theta), when
% accounting for the gradient of P(o=1|k',theta) wrt theta. But if the
% agent is uncertain about k', then she will rescale the update according
% to P(k'), ie the effective prediction error becomes P(k')*PE.
% Note that k-ToM assumes that her opponent's parameters can drift over
% trials. How much they can change is controlled by her prior volatility,
% which is the only unknown parameter of k-ToM's evolution function.

a = 0.36; % for E[s(x)] when x~n(mu,Sig)
level = inF.lev; % depth of k-ToM's recursive beliefs
fx = zeros(size(x)); % initialize hidden states

if level==0 % 0-ToM
    
    fx = evolution0bisND(x,P,u,inF);  % here: simple Laplace-Kalman filter
    indlev = [];
    
else % k-ToM with k>0
    
    fobs = inF.fobs; % opponent's observation function (ie decision mapping)
    NParev = inF.indParev; % number of opponent's evolution params
    NParobs = inF.indParobs;  % number of opponent's observation params
    NtotPar = inF.indParev+inF.indParobs; % total number of params
    indlev = defIndlev(level,NtotPar); % states indexing
    
    % 1- Update P(k'), the agent's belief Re: her opponent's level
    % Recall that k-ToM maintains parallel predictions P(o=1|k')
    % about her opponent's next move, conditional upon her opponent's level
    % k'. In brief, if P(o=1|k') ~= o, then this is evidence that her
    % opponent has level k'. In turn, P(k') is increased if P(o=1|k')
    % is the closest to o (across different k')...
    % Note: hidden states encode x(theta), the log-odds of P(o=1|k',theta)
    % evaluated at the agent's estimate of theta (mu). In addition, they
    % encode the gradient of x wrt to theta (dx/dtheta), and V[theta]. This
    % then serves to derive P(o=1|k') as follows:
    % P(o=1|k') = E[sigm(x(theta))]
    %           = sigm(E[x(theta)]/sqrt(1+a*V[x(theta)])
    %           = sigm(E[x(theta)]/sqrt(1+a*V[theta]*(dx/dtheta)^2)
    if isempty(u) % missed trial or initialization
        % max-entropic belief about opponent's level
        Pu = 1./level.*ones(1,level);
        if level>1
            fx(1:(level-1)) = invsigmoid(Pu(1:(end-1)));
        end
    else % trial OK
        y2t = u(1); % opponent's last move
        if level<2 % 1-ToM
            Pu = 1; % 1-ToM's opponent = 0-ToM!
        else % k-ToM with k>=2
            % get sufficient statistics of posterior belief on x(theta)
            f = zeros(1,level); % E[x(theta)]
            Vx = zeros(1,level); % V[x(theta)]
            for j=1:level % loop over possible opponent's levels  (k'=j-1)
                f(j) = x(indlev(j).f); % E[x(theta)|k'=j-1]
                df = x(indlev(j).df); % d[x(theta)]/dtheta for k'=j-1
                Sig = exp(x(indlev(j).Par(2:2:2*NtotPar))); % V[theta|k'=j-1]
                Vx(j) = Sig'*df.^2;
            end
            % derive E[sigm(x(theta))]
            Es = sigmoid(f./sqrt(1+a*Vx)); % ~= MD's derivation!!!
            % get prior P(k')
            P0 = zeros(1,level);
            P0(1:(level-1)) = sigmoid(x(1:(level-1)));
            P0(end) = max(.0001,1-sum(P0)); % insert last P(k=k'-1)
            % partial forgetting of opponent's level?
            if inF.diluteP % this is not in [Devaine et al. (20014)]
                dc = sigmoid(P(2)); % dilution coefficient
                P0 = (1-dc).*P0 + dc./(level.*ones(1,level));
            end
            % update P(k')
            Pu = P0.*Es.^y2t.*(1-Es).^(1-y2t);
            Pu = Pu./sum(Pu); % posterior P(k')
            % store posterior P(k') in hidden-states vector
            fx(1:(level-1)) = invsigmoid(Pu(1:(end-1))); % only k'<k-1
        end
    end
    
    
    % 2- Update P(theta|k'), the agent's belief Re: her opponent's params
    % Note: the opponent's parameters can either be assumed to drift over
    % trials, or be (a priori) invariant. In case a parameter is assumed to
    % drift, its previous posterior belief is diluted in proportion to
    % exp(P), which is k-ToM's prior volatility...
    for j=1:level % loop over admissible opponent's levels (k'=j-1)
        
        % 2.1- update E[theta|k'] and V[theta|k']
        in.k = Pu(j); % posterior P(k')
        in.f = x(indlev(j).f); % E[x(theta)|k'=j]
        in.df = x(indlev(j).df); % d[x(theta)]/dtheta for k'=j
        in.dummyPar = inF.dummyPar; % flags which params drift over trials
        if isempty(u) % missed trial or initialization
            fx(indlev(j).Par) = x(indlev(j).Par); % no update
            u_t = [];
        else
            fx(indlev(j).Par) = updatePar(x(indlev(j).Par),P(1),y2t,in);
            u_t = u([2;1]); % reverse inputs 
        end
        %-----------------------------------------------------------------%
        % Note: this concludes k-ToM's update rule. What follows is not
        % formally mandatory, but is justified from practical reasons. In
        % brief, k-ToM's hidden states are augmented with dummy states that
        % are derived from sufficent statistics of posterior beliefs P(k')
        % and p(theta).
        %-----------------------------------------------------------------%
        
        inL = inF; % copy optional input structure
        inL.lev = j-1; % opponent sophistication [j=1:0-ToM, j=2:1-ToM,...]
        inL.player = 3-inF.player; % reverse player role
        
        % 2.2- simulate evolution of opponent's hidden states 
        indMu_ev = indlev(j).Par(1:2:2*NParev); % indices of opponent's evol params
        indMu_obs = indlev(j).Par(2*NParev+1:2:2*NtotPar); % indices of opponent's obs params
        Parjev = fx(indMu_ev); % E[evol param 
        Parobs = fx(indMu_obs); % E[obs param]
        [X,indlevj] = RecToMfunction(x(indlev(j).X),Parjev,u_t,inL);
        fx(indlev(j).X) = X; % store evolving states
        
        % 2.3- update log-odds of P(o=1|k'=j-1), ie x(params)
        inG = inF; % copy optional input structure
        inG.npara = NtotPar; % total nb of evol/obs params
        inG.lev = j-1; % opponent sophistication
        inG.k = in.k; % [useless]
        inG.player = 3-inF.player; % reverse player role
        inG.indlev = indlevj; % indexing of opponent's hidden states
        f_new = invsigmoid(fobs(X,Parobs,u_t,inG)); % x(params)
        fx(indlev(j).f) = f_new; % store evolving states
        
        % 2.4- get numerical derivative of x(params) wrt evol params
        eps = zeros(1,NParev);
        f_neweps = zeros(NParev,1);
        df_newev = zeros(NParev,1);
        for l=1:NParev % loop over evol params
            eps(l) = 1e-4*Parjev(l); % small param increment (dP)
            if abs(eps(l))<1e-4
                eps(l) = 1e-4;
            end
            Parjevl = Parjev;
            Parjevl(l) = Parjevl(l)+eps(l); % P <-- P + dP
            Xeps = RecToMfunction(x(indlev(j).X),Parjevl,u_t,inL);
            f_neweps(l) = invsigmoid(fobs(Xeps,Parobs,u_t,inG));
            df_newev(l) = (f_neweps(l)-f_new)/eps(l);
        end
        fx(indlev(j).df(1:NParev)) = df_newev;
        
        % 2.5- get numerical derivative of x(params) wrt obs params
        eps = zeros(1,NParobs);
        f_neweps = zeros(NParobs,1);
        df_newobs = zeros(NParobs,1);
        for l=1:NParobs % loop over obs params
            eps(l)=1e-4*Parobs(l); % small param increment (dP)
            if abs(eps(l))<1e-4
                eps(l) = 1e-4;
            end
            Parobsl = Parobs;
            Parobsl(l) = Parobsl(l)+eps(l); % P <-- P + dP
            f_neweps(l) = invsigmoid(fobs(X,Parobsl,u_t,inG));
            df_newobs(l) = (f_neweps(l)-f_new)/eps(l);
        end
        fx(indlev(j).df(NParev+(1:NParobs))) = df_newobs;
        
    end
end
end % for avoiding shared variables in nested functions 


function fPar_k = updatePar(Par_k,theta,u_t,in)
% NB: VB update of opponent's params does not account for potential
% identifiability issues between params (cf. mean-field assumption)
Pu = in.k; % updated P(k)
Proba = sigmoid(in.f); % P(o=1|k)
df = in.df; % d[x(theta)]/dtheta
indV = 2:2:length(Par_k); % indices of E[theta]
indMu = 1:2:length(Par_k); % indices of V[theta]
V0 = exp(Par_k(indV))+exp(theta).*in.dummyPar; % diluted prior variance
Vu = 1./((1./V0)+Pu*Proba*(1-Proba).*(df.^2)); % posterior variance
E0 = Par_k(indMu); % prior mean
Eu = E0 + Pu.*Vu.*(u_t-Proba).*df; % posterior mean
fPar_k = zeros(size(Par_k)); % = [...,E[param_i];V[param_i];...]
fPar_k(indV) = log(Vu);
fPar_k(indMu) = Eu;
end
